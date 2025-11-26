#pragma once
#include <vector>
#include <string>
#include <filesystem>


struct Triangle {
    std::string ID;
    double v1_x, v1_y, v1_z; // Вершина 1
    double v2_x, v2_y, v2_z; // Вершина 2
    double v3_x, v3_y, v3_z; // Вершина 3
    double normal_x, normal_y, normal_z; // Нормаль
    std::string component_type; // Тип компонента (Body, Panel, Antenna)
    double label; // Метка (1 - освещен, 0 - в тени)
};

class SatelliteDataset {
public:
    SatelliteDataset(const std::string& folder_path);
    const std::vector<Triangle>& get_triangles() const { return triangles; }
    void load_single_file(const std::string& file_path);
    const std::vector<std::string>& get_file_names() const { return file_names; }

private:
    std::vector<Triangle> triangles;
    std::vector<std::string> file_names;

    void load_data(const std::string& folder_path);
};