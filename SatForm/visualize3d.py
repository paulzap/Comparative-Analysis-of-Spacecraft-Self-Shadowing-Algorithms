import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# === Цвета ===
colors = {
    'yellow': (250/255, 203/255, 105/255),
    'red':    (222/255, 95/255,  100/255),
    'green':  (100/255, 175/255, 135/255),
    'orange': (250/255, 157/255, 91/255),
    'blue':   (95/255,  183/255, 230/255),
    'violet': (165/255, 140/255, 190/255),
}

component_color = {
    "Body":    colors['blue'],
    "Panel":   colors['green'],
    "Antenna": colors['red'],
}

def visualize_triangles(
    csv_file,
    step=1,
    show_normals=False,
    alpha=0.8,
    linewidth=0.07,
    normal_scale=0.8
):
    logging.info(f"Чтение файла: {csv_file}")
    df = pd.read_csv(csv_file)

    # Принудительно переводим нужные столбцы в float (это решает все проблемы с типами)
    coord_cols = [f"V{i}_{c}" for i in (1,2,3) for c in ("X","Y","Z")]
    normal_cols = ["Normal_X", "Normal_Y", "Normal_Z"]
    all_numeric_cols = coord_cols + normal_cols + ["Triangle ID"]  # ID тоже может мешать

    for col in all_numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')  # всё что не число → NaN

    # Удаляем строки, где есть NaN в координатах
    df = df.dropna(subset=coord_cols)

    if df.empty:
        logging.error("После очистки данных не осталось ни одной валидной строки!")
        return

    verts_list   = []
    facecolors   = []
    centers      = []
    normals_vec  = []
    panel_hashes = set()           # только для панелей — убираем обратные грани

    for comp_type in ["Body", "Panel", "Antenna"]:
        subset = df[df["Component Type"] == comp_type].iloc[::step]
        if subset.empty:
            continue

        # Берём только числовые столбцы в правильном порядке
        verts_arr   = subset[coord_cols].values.reshape(-1, 3, 3)   # [N,3,3]
        normals_arr = subset[normal_cols].values                    # [N,3]

        added = 0
        for i in range(len(verts_arr)):
            tri = verts_arr[i]          # shape (3,3)
            nrm = normals_arr[i]        # shape (3,)

            # Пропускаем полностью битые треугольники
            if np.any(np.isnan(tri)) or np.any(np.isnan(nrm)):
                continue

            # === Убираем обратные грани только у панелей ===
            if comp_type == "Panel":
                # Канонический хэш: сортируем 3 вершины по координатам
                rounded = np.round(tri, decimals=6)
                # Превращаем в список кортежей и сортируем
                verts_tuple = tuple(tuple(rounded[j]) for j in range(3))
                sorted_tuple = tuple(sorted(verts_tuple))
                h = hash(sorted_tuple)

                if h in panel_hashes:
                    continue
                panel_hashes.add(h)

            # Добавляем в итоговый список
            verts_list.append(tri)
            facecolors.append(component_color[comp_type])
            added += 1

            if show_normals:
                centers.append(tri.mean(axis=0))
                normals_vec.append(nrm)

        logging.info(f"Добавлено {added} треугольников типа {comp_type}")

    if not verts_list:
        logging.error("После фильтрации ничего не осталось!")
        return

    # === Рисуем ===
    fig = plt.figure(figsize=(11, 10))
    ax = fig.add_subplot(111, projection='3d')

    collection = Poly3DCollection(
        verts_list,
        facecolors=facecolors,
        edgecolors='black',
        linewidths=linewidth,
        alpha=alpha,
        shade=True,
        antialiased=True,
    )
    ax.add_collection3d(collection)

    # Нормали
    if show_normals and centers:
        centers = np.array(centers)
        normals_vec = np.array(normals_vec)
        lengths = np.linalg.norm(normals_vec, axis=1, keepdims=True)
        lengths[lengths == 0] = 1.0
        normals_vec = normals_vec / lengths * normal_scale

        ax.quiver(
            centers[:,0], centers[:,1], centers[:,2],
            normals_vec[:,0], normals_vec[:,1], normals_vec[:,2],
            length=1.0, normalize=True,
            color='white', alpha=0.8, linewidth=1.2
        )

    # === Кубический масштаб ===
    all_v = np.concatenate(verts_list, axis=0)
    mins, maxs = all_v.min(axis=0), all_v.max(axis=0)
    center = (mins + maxs) / 2
    r = (maxs - mins).max() / 2 * 1.05

    ax.set_xlim(center[0] - r, center[0] + r)
    ax.set_ylim(center[1] - r, center[1] + r)
    ax.set_zlim(center[2] - r, center[2] + r)
    ax.set_box_aspect([1, 1, 1])

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Спутник — односторонние панели, корректная глубина')

    ax.xaxis.pane.fill = ax.yaxis.pane.fill = ax.zaxis.pane.fill = False
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()


# === Запуск ===
# visualize_triangles("your_triangles.csv", step=1, show_normals=False, alpha=0.9, linewidth=0.2)