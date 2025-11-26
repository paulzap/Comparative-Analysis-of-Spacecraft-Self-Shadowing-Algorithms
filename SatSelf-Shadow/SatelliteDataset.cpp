#include "SatelliteDataset.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

SatelliteDataset::SatelliteDataset(const std::string& folder_path) {
    load_data(folder_path);
}

void SatelliteDataset::load_data(const std::string& folder_path) {
    // Получаем список CSV-файлов
    file_names.clear();
    for (const auto& entry : std::filesystem::directory_iterator(folder_path)) {
        if (entry.path().extension() == ".csv") {
            file_names.push_back(entry.path().string());
        }
    }
    if (file_names.empty()) {
        throw std::runtime_error("No CSV files found in " + folder_path);
    }

    // Сортируем файлы для воспроизводимости
    std::sort(file_names.begin(), file_names.end());
}

void SatelliteDataset::load_single_file(const std::string& file_path) {
    triangles.clear();
    std::ifstream file(file_path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open " + file_path);
    }

    // Пропускаем заголовок
    std::string line;
    std::getline(file, line); // Пропускаем "Triangle ID,Component Type,..."

    // Читаем данные
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        Triangle triangle;
        triangle.label = 1.0; // Инициализация метки

        // Читаем значения, разделенные запятыми
        std::getline(ss, token, ','); triangle.ID = token;// Пропускаем Triangle ID
        std::getline(ss, token, ','); triangle.component_type = token;
        std::getline(ss, token, ','); triangle.v1_x = std::stod(token);
        std::getline(ss, token, ','); triangle.v1_y = std::stod(token);
        std::getline(ss, token, ','); triangle.v1_z = std::stod(token);
        std::getline(ss, token, ','); triangle.v2_x = std::stod(token);
        std::getline(ss, token, ','); triangle.v2_y = std::stod(token);
        std::getline(ss, token, ','); triangle.v2_z = std::stod(token);
        std::getline(ss, token, ','); triangle.v3_x = std::stod(token);
        std::getline(ss, token, ','); triangle.v3_y = std::stod(token);
        std::getline(ss, token, ','); triangle.v3_z = std::stod(token);
        std::getline(ss, token, ','); triangle.normal_x = std::stod(token);
        std::getline(ss, token, ','); triangle.normal_y = std::stod(token);
        std::getline(ss, token, ','); triangle.normal_z = std::stod(token);

        // Нормализация нормали
        double length = std::sqrt(triangle.normal_x * triangle.normal_x +
            triangle.normal_y * triangle.normal_y +
            triangle.normal_z * triangle.normal_z);
        if (length > 0.0) { // Избегаем деления на ноль
            triangle.normal_x /= length;
            triangle.normal_y /= length;
            triangle.normal_z /= length;
        }

        triangles.push_back(triangle);
    }

    file.close();
}