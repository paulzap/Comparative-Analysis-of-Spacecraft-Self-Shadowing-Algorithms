#pragma once
#include "CLASS_SpacecraftPart.h"
#include <memory>
#include <vector>
#include <fstream>
#include <iomanip>

class Spacecraft {
private:
    std::unique_ptr<SpacecraftPart> body;
    std::vector<std::unique_ptr<SpacecraftPart>> solarPanels;
    std::unique_ptr<SpacecraftPart> antenna;

public:
    Spacecraft(double _pol_size) {
        body = std::make_unique<Body>();
        body->initialize(_pol_size);

        body->generateVerticies();

        std::vector<Point> limits = body->getShape().getExtremePoints();

        // First solar panel (always present)
        auto panel1 = std::make_unique<SolarPanel>();
        panel1->initialize(_pol_size);
        panel1->generateVerticies(limits, 0);
        solarPanels.push_back(std::move(panel1));

        // Second solar panel (75% probability)
        int randValue = rand() % 100;
        if (randValue >= 25) {
            auto panel2 = std::make_unique<SolarPanel>();
            panel2->initialize(_pol_size);
            panel2->generateVerticies(limits, 2);
            solarPanels.push_back(std::move(panel2));
        }

        antenna = std::make_unique<Antenna>();
        antenna->initialize(_pol_size);
        antenna->generateVerticies(limits, 3);
    }

    void saveToFile(std::ofstream& out, int satelliteNumber = 1) const {
        out << "Satellite #" << satelliteNumber << "\n";
        body->storeToFile(out);
        for (const auto& panel : solarPanels) {
            panel->storeToFile(out);
        }
        if (antenna) {
            antenna->storeToFile(out);
        }
        out << "-------------------------\n";
    }

    std::tuple<std::vector<Triangle>, std::vector<Triangle>, std::vector<Triangle>> getTriangles() const {
        std::vector<Triangle> body_triangles = body->getTriangles();
        std::vector<Triangle> panel_triangles;
        for (const auto& panel : solarPanels) {
            auto panel_tris = panel->getTriangles();
            panel_triangles.insert(panel_triangles.end(), panel_tris.begin(), panel_tris.end());
        }
        std::vector<Triangle> antenna_triangles = antenna ? antenna->getTriangles() : std::vector<Triangle>();
        return { body_triangles, panel_triangles, antenna_triangles };
    }

    void partitionShapes( const std::string& filename) {
        auto [body_triangles, panel_triangles, antenna_triangles] = getTriangles();
        saveTrianglesToCSV(body_triangles, panel_triangles, antenna_triangles, filename);
    }

    void saveTrianglesToCSV(const std::vector<Triangle>& body_triangles, const std::vector<Triangle>& panel_triangles,
        const std::vector<Triangle>& antenna_triangles, const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file for writing: " + filename);
        }

        // Устанавливаем максимальную точность для double
        file << std::setprecision(std::numeric_limits<double>::max_digits10);

        // Заголовок CSV
        file << "Triangle ID,Component Type,V1_X,V1_Y,V1_Z,V2_X,V2_Y,V2_Z,V3_X,V3_Y,V3_Z,Normal_X,Normal_Y,Normal_Z\n";
        int id = 0;

        // Функция для форматирования чисел: замена близких к нулю значений и устранение -0
        auto format_number = [](double value) -> double {
            const double epsilon = 1e-15; // Порог для замены на 0
            if (std::abs(value) < epsilon) {
                return 0.0; // Заменяем -0 и e-17, e-18 на 0
            }
            return value;
            };

        // Запись треугольников тела
        for (const auto& triangle : body_triangles) {
            auto normals = body->getShape().computeTriangleNormals(triangle);
            for (const auto& normal : normals) {
                file << id++ << ",Body,"
                    << format_number(triangle.v1.x) << "," << format_number(triangle.v1.y) << "," << format_number(triangle.v1.z) << ","
                    << format_number(triangle.v2.x) << "," << format_number(triangle.v2.y) << "," << format_number(triangle.v2.z) << ","
                    << format_number(triangle.v3.x) << "," << format_number(triangle.v3.y) << "," << format_number(triangle.v3.z) << ","
                    << format_number(normal.x) << "," << format_number(normal.y) << "," << format_number(normal.z) << "\n";
            }
        }

        // Запись треугольников панелей
        for (const auto& triangle : panel_triangles) {
            auto normals = solarPanels[0]->getShape().computeTriangleNormals(triangle);
            for (const auto& normal : normals) {
                file << id++ << ",Panel,"
                    << format_number(triangle.v1.x) << "," << format_number(triangle.v1.y) << "," << format_number(triangle.v1.z) << ","
                    << format_number(triangle.v2.x) << "," << format_number(triangle.v2.y) << "," << format_number(triangle.v2.z) << ","
                    << format_number(triangle.v3.x) << "," << format_number(triangle.v3.y) << "," << format_number(triangle.v3.z) << ","
                    << format_number(normal.x) << "," << format_number(normal.y) << "," << format_number(normal.z) << "\n";
            }
        }

        // Запись треугольников антенны
        for (const auto& triangle : antenna_triangles) {
            auto normals = antenna->getShape().computeTriangleNormals(triangle);
            for (const auto& normal : normals) {
                file << id++ << ",Antenna,"
                    << format_number(triangle.v1.x) << "," << format_number(triangle.v1.y) << "," << format_number(triangle.v1.z) << ","
                    << format_number(triangle.v2.x) << "," << format_number(triangle.v2.y) << "," << format_number(triangle.v2.z) << ","
                    << format_number(triangle.v3.x) << "," << format_number(triangle.v3.y) << "," << format_number(triangle.v3.z) << ","
                    << format_number(normal.x) << "," << format_number(normal.y) << "," << format_number(normal.z) << "\n";
            }
        }

        file.close();
    }
};
