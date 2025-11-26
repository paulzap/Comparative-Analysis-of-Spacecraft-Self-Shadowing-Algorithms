import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import Axes3D
# Явное указание бэкенда
try:
    matplotlib.use('Qt5Agg')
except ImportError:
    try:
        matplotlib.use('TkAgg')
    except ImportError:
        pass

def visualize_triangles(csv_file, step=1, show_normals=False):
    """
    Визуализирует треугольники спутника из CSV-файла в 3D с использованием Matplotlib.
    Треугольники окрашиваются в зависимости от метки освещенности: желтый для освещенных (Label=1),
    серый для затененных (Label=0). Клоны треугольников (с одинаковыми вершинами и противоположными
    нормалями) сдвигаются вдоль нормали для устранения z-fighting, оригинал скрывается.
    Границы треугольников не отображаются для уменьшения визуального шума. Шаг выборки автоматически
    корректируется для больших сцен.
    Parameters:
    -----------
    csv_file : str
        Путь к CSV-файлу с данными треугольников.
    step : int
        Шаг выборки треугольников (например, 1 - все треугольники, 10 - каждый 10-й).
        Если step=0, используется автоматический шаг для ограничения числа полигонов.
    show_normals : bool
        Если True, отображает нормали для треугольников (по умолчанию False).
    """
    # Загрузка данных из CSV
    try:
        df = pd.read_csv(csv_file)
    except FileNotFoundError:
        print(f"File {csv_file} not found")
        return
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return
    # Извлечение треугольников, нормалей и меток
    triangles = df[["V1_X", "V1_Y", "V1_Z", "V2_X", "V2_Y", "V2_Z", "V3_X", "V3_Y", "V3_Z"]].values.reshape(-1, 3, 3)
    normals = df[["Normal_X", "Normal_Y", "Normal_Z"]].values
    labels = df["Label"].values
    correct = df["Correct"].values if "Correct" in df.columns else np.ones(len(labels), dtype=int)
    if len(triangles) == 0:
        print("No triangles found in CSV file, nothing to visualize")
        return
    # Автоматический выбор шага
    max_polygons = 5000 # Максимальное число полигонов для рендеринга
    if step == 0:
        step = max(1, len(triangles) // max_polygons + 1)
        if step > 1:
            print(f"Using automatic step={step} to limit polygons to ~{max_polygons}")
    # Обнаружение клонов с хэшированием и формирование видимых треугольников
    epsilon = 0.01  # Расстояние сдвига (увеличено для видимости)
    tolerance = 1e-6  # Допуск для сравнения
    vertex_hash = {}  # Хэш-таблица: ключ - вершины (округленные), значение - список индексов
    for i in range(len(triangles)):
        # Округление вершин для хэширования
        verts = tuple(np.round(triangles[i].flatten() / tolerance).astype(int))
        if verts not in vertex_hash:
            vertex_hash[verts] = []
        vertex_hash[verts].append(i)
    
    # Списки для видимых треугольников и меток
    visible_triangles = []
    visible_labels = []
    visible_normals = []  # Для нормалей, если нужно
    visible_correct = []  # For correctness tracking
    for indices in vertex_hash.values():
        if len(indices) == 1:
            # Одиночный треугольник: добавляем как есть
            i = indices[0]
            visible_triangles.append(triangles[i])
            visible_labels.append(labels[i])
            visible_correct.append(correct[i])
            visible_normals.append(normals[i])
        else:
            # Группа с клонами: ищем пары с противоположными нормалями
            processed = set()
            for i in indices:
                if i in processed:
                    continue
                for j in indices:
                    if i == j or j in processed:
                        continue
                    # Проверка противоположных нормалей
                    normal_opposite = np.all(np.abs(normals[i] + normals[j]) < tolerance)
                    if normal_opposite:
                        # Пара найдена: добавляем два сдвинутых треугольника, оригинал скрываем
                        verts_base = triangles[i]  # Базовые вершины (одинаковые)
                        # Сдвинутый от первого (оригинал)
                        shifted1 = verts_base + epsilon * normals[i]
                        visible_triangles.append(shifted1)
                        visible_labels.append(labels[i])
                        visible_correct.append(correct[i])
                        visible_normals.append(normals[i])
                        # Сдвинутый от второго (клон)
                        shifted2 = verts_base + epsilon * normals[j]
                        visible_triangles.append(shifted2)
                        visible_labels.append(labels[j])
                        visible_correct.append(correct[j])
                        visible_normals.append(normals[j])
                        processed.add(i)
                        processed.add(j)
                        break  # Переходим к следующей паре
                if i not in processed:
                    # Если нет пары, добавляем как одиночный (сдвиг, если нужно)
                    visible_triangles.append(triangles[i])
                    visible_labels.append(labels[i])
                    visible_correct.append(correct[i])
                    visible_normals.append(normals[i])
                    processed.add(i)
    
    # Преобразование в массивы
    visible_triangles = np.array(visible_triangles)
    visible_labels = np.array(visible_labels)
    visible_normals = np.array(visible_normals)
    visible_correct = np.array(visible_correct)
    # Применение шага выборки
    step_indices = list(range(0, len(visible_triangles), step))
    verts = visible_triangles[step_indices]
    labels_step = visible_labels[step_indices]
    correct_step = visible_correct[step_indices]
    colors = ['red' if correct_step[i] == 0 else ('yellow' if labels_step[i] == 1 else 'gray') for i in range(len(labels_step))]
    
    # Создание фигуры и 3D-осей
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    # Визуализация треугольников
    collection = Poly3DCollection(verts, alpha=1.0, facecolors=colors, edgecolor='none')
    ax.add_collection3d(collection)
    # Визуализация нормалей (если включено)
    if show_normals and len(visible_normals) > 0:
        normals_step = visible_normals[step_indices]
        centers = np.mean(verts, axis=1)
        normal_colors = ['yellow' if label == 1 else 'gray' for label in labels_step]
        for center, normal, color in zip(centers, normals_step, normal_colors):
            ax.quiver(center[0], center[1], center[2],
                      normal[0], normal[1], normal[2],
                      color=color, length=1.0, normalize=True)
    # Настройка одинакового масштаба осей (используем все вершины, включая сдвинутые)
    all_verts = np.vstack([visible_triangles, triangles])  # Для полного диапазона
    min_coords = all_verts.min(axis=(0, 1))
    max_coords = all_verts.max(axis=(0, 1))
    max_range = np.max(max_coords - min_coords)
    mid_coords = (max_coords + min_coords) / 2
    ax.set_xlim(mid_coords[0] - max_range / 2, mid_coords[0] + max_range / 2)
    ax.set_ylim(mid_coords[1] - max_range / 2, mid_coords[1] + max_range / 2)
    ax.set_zlim(mid_coords[2] - max_range / 2, mid_coords[2] + max_range / 2)
    ax.set_box_aspect([1, 1, 1])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('3D Visualization of Satellite with Shadows (Originals Hidden)')
    # Отображение графика
    plt.show(block=True)
