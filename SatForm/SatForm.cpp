#include "stdafx.h"
#include "CLASS_Spacecraft.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <Python.h>

namespace py = pybind11;

static struct PyInitializer {
    PyInitializer() {
        PyConfig config;
        PyConfig_InitPythonConfig(&config);
        PyStatus status = PyConfig_SetString(&config, &config.program_name, L"C:\\Users\\zapevalin\\AppData\\Local\\Programs\\Python\\Python312\\python.exe");
        if (PyStatus_IsError(status)) {
            throw std::runtime_error("Failed to set program_name: " + std::string(status.err_msg));
        }
        status = PyConfig_SetString(&config, &config.executable, L"C:\\Users\\zapevalin\\AppData\\Local\\Programs\\Python\\Python312\\python.exe");
        if (PyStatus_IsError(status)) {
            throw std::runtime_error("Failed to set executable: " + std::string(status.err_msg));
        }
        status = PyConfig_SetString(&config, &config.home, L"C:\\Users\\zapevalin\\AppData\\Local\\Programs\\Python\\Python312");
        if (PyStatus_IsError(status)) {
            throw std::runtime_error("Failed to set PYTHONHOME: " + std::string(status.err_msg));
        }
        config.module_search_paths_set = 1;
        status = PyWideStringList_Append(&config.module_search_paths, L"C:\\Users\\zapevalin\\AppData\\Local\\Programs\\Python\\Python312\\Lib");
        if (PyStatus_IsError(status)) {
            throw std::runtime_error("Failed to append Python Lib path: " + std::string(status.err_msg));
        }
        status = PyWideStringList_Append(&config.module_search_paths, L"C:\\Users\\zapevalin\\AppData\\Local\\Programs\\Python\\Python312\\Lib\\site-packages");
        if (PyStatus_IsError(status)) {
            throw std::runtime_error("Failed to append site-packages path: " + std::string(status.err_msg));
        }
        status = PyWideStringList_Append(&config.module_search_paths, L"C:\\Users\\zapevalin\\AppData\\Local\\Programs\\Python\\Python312\\DLLs");
        if (PyStatus_IsError(status)) {
            throw std::runtime_error("Failed to append DLLs path: " + std::string(status.err_msg));
        }
        status = PyWideStringList_Append(&config.module_search_paths, L"e:\\Ïåðñïåêòèâíûå çàäà÷è\\Çàäà÷à 18. Íåéðîñåòü òåíåé\\Ìîäåëèðîâàíèå\\SatForm3D\\SatForm");
        if (PyStatus_IsError(status)) {
            throw std::runtime_error("Failed to append project path: " + std::string(status.err_msg));
        }
        status = Py_InitializeFromConfig(&config);
        if (PyStatus_IsError(status)) {
            throw std::runtime_error("Failed to initialize Python: " + std::string(status.err_msg));
        }
        PyConfig_Clear(&config);
    }
    ~PyInitializer() { py::finalize_interpreter(); }
} py_init;

void visualize_triangles(const std::string& csv_file, int step = 1) {
    try {
        py::module_ visualize = py::module_::import("visualize3d");
        visualize.attr("visualize_triangles")(csv_file, step);
    }
    catch (const py::error_already_set& e) {
        std::cerr << "Python error: " << e.what() << "\n";
        throw;
    }
}

#include <filesystem> // Äîáàâüòå ýòîò include

void generateSatellites(int count, double _pol_size, const std::string& folder, const std::string& filename, int step = 1) {
    // Ñîçäàåì ïàïêó, åñëè îíà íå ñóùåñòâóåò
    if (!std::filesystem::exists(folder)) {
        std::filesystem::create_directories(folder);
        std::cout << "Created directory: " << folder << "\n";
    }

    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Failed to open file: " << filename << "\n";
        return;
    }

    for (int i = 1; i <= count; ++i) {
        Spacecraft sat(_pol_size);
        sat.saveToFile(out, i);
        std::string csv_filename = folder + "/spc_" + std::to_string(i) + "_" + std::to_string(_pol_size).substr(0, 3) + ".csv";
        sat.partitionShapes(csv_filename);

        // Âèçóàëèçàöèÿ íà îñíîâå CSV-ôàéëà
        //if (i == -1)
            visualize_triangles(csv_filename, step);
        //}

        if (i % 100 == 0 || i == count) {
            std::cout << std::setw(6) << i << "/" << count << " complete (" << static_cast<double>(i) / count * 100 << "%)" << std::endl;
        }
    }

    out.close();
}

int main() {
    // std::srand(static_cast<unsigned int>(std::time(nullptr)));
    //generateSatellites(1000, 1.5, "satellites.txt", 1); // pol_size=0.5, step=2
    
    
    int num = 1000;
    std::string folder = "data3d_1.5";
    generateSatellites(num, 1.5, folder, folder + "/satellites.txt", 1);

    //folder = "data3d_0.05";
    //generateSatellites(num, 0.05, folder, folder + "/satellites.txt", 1);

    //folder = "data3d_0.2";
    //generateSatellites(num, 0.2, folder, folder + "/satellites.txt", 1);

    return 0;

}
