#include "SatelliteDataset.h"
#include "ShadowAlgorithms.h"
#include <iostream>
#include <chrono>
#include <vector>
#include <string>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <random>
#include <cmath>
#include <fstream>
#include <filesystem>

#define M_PI 3.14159265358979323846 /* pi */

namespace fs = std::filesystem;

// Structure to store benchmark results for each run
struct BenchmarkResult {
    std::string method_name;
    std::string file_name;
    int run_number;
    double time_ms;
    double accuracy;
    double dice;
    double jaccard;
    double theta_degrees;
    double phi_degrees;
    double alpha_degrees;
    int num_triangles;
    int true_positive;
    int true_negative;
    int false_positive;
    int false_negative;
    double f1_score;
    double precision;
    double recall;
};

// Structure to hold confusion matrix metrics
struct ConfusionMetrics {
    int true_positive = 0;
    int true_negative = 0;
    int false_positive = 0;
    int false_negative = 0;
    double f1_score = 0.0;
    double precision = 0.0;
    double recall = 0.0;
};

// Calculate confusion matrix metrics
ConfusionMetrics calculate_confusion_metrics(const std::vector<int>& predicted, const std::vector<int>& reference) {
    ConfusionMetrics metrics;
    if (predicted.size() != reference.size() || predicted.empty()) return metrics;
    for (size_t i = 0; i < predicted.size(); ++i) {
        if (predicted[i] == 1 && reference[i] == 1) metrics.true_positive++;
        else if (predicted[i] == 0 && reference[i] == 0) metrics.true_negative++;
        else if (predicted[i] == 1 && reference[i] == 0) metrics.false_positive++;
        else if (predicted[i] == 0 && reference[i] == 1) metrics.false_negative++;
    }
    if (metrics.true_positive + metrics.false_positive > 0)
        metrics.precision = static_cast<double>(metrics.true_positive) / (metrics.true_positive + metrics.false_positive);
    if (metrics.true_positive + metrics.false_negative > 0)
        metrics.recall = static_cast<double>(metrics.true_positive) / (metrics.true_positive + metrics.false_negative);
    if (metrics.precision + metrics.recall > 0)
        metrics.f1_score = 2.0 * (metrics.precision * metrics.recall) / (metrics.precision + metrics.recall);
    return metrics;
}



// Calculate accuracy (percentage of correctly classified labels)
double calculate_accuracy(const std::vector<int>& predicted, const std::vector<int>& reference) {
    if (predicted.size() != reference.size() || predicted.empty()) {
        return 0.0;
    }
    size_t matches = 0;
    for (size_t i = 0; i < predicted.size(); ++i) {
        if (predicted[i] == reference[i]) {
            matches++;
        }
    }
    return (matches / static_cast<double>(predicted.size())) * 100.0;
}

// Calculate Dice coefficient (F1 score)
double calculate_dice(const std::vector<int>& predicted, const std::vector<int>& reference) {
    if (predicted.size() != reference.size() || predicted.empty()) {
        return 0.0;
    }
    
    int true_positive = 0;
    int false_positive = 0;
    int false_negative = 0;
    
    for (size_t i = 0; i < predicted.size(); ++i) {
        if (predicted[i] == 1 && reference[i] == 1) {
            true_positive++;
        } else if (predicted[i] == 1 && reference[i] == 0) {
            false_positive++;
        } else if (predicted[i] == 0 && reference[i] == 1) {
            false_negative++;
        }
    }
    
    if (true_positive + false_positive + false_negative == 0) {
        return 100.0; // If all are true negatives
    }
    
    return ((2.0 * true_positive) / (2.0 * true_positive + false_positive + false_negative)) * 100.0;
}

// Calculate Jaccard index (Intersection over Union)
double calculate_jaccard(const std::vector<int>& predicted, const std::vector<int>& reference) {
    if (predicted.size() != reference.size() || predicted.empty()) {
        return 0.0;
    }
    
    int intersection = 0;
    int union_count = 0;
    
    for (size_t i = 0; i < predicted.size(); ++i) {
        if (predicted[i] == 1 && reference[i] == 1) {
            intersection++;
        }
        if (predicted[i] == 1 || reference[i] == 1) {
            union_count++;
        }
    }
    
    if (union_count == 0) {
        return 100.0; // If all are zeros (both predicted and reference)
    }
    
    return (static_cast<double>(intersection) / union_count) * 100.0;
}

// Save benchmark results to CSV
void save_benchmark_results(const std::vector<BenchmarkResult>& results, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }
    
    // Write header
    file << "method_name,file_name,run_number,num_triangles,time_ms,accuracy,dice,jaccard,theta_degrees,phi_degrees,alpha_degrees,true_positive,true_negative,false_positive,false_negative,f1_score,precision,recall\n";
    
    // Write data
    for (const auto& result : results) {
        file << result.method_name << ","
             << result.file_name << ","
             << result.run_number << ","
             << result.num_triangles << ","
             << std::fixed << std::setprecision(3) << result.time_ms << ","
             << std::fixed << std::setprecision(6) << result.accuracy << ","
             << std::fixed << std::setprecision(6) << result.dice << ","
             << std::fixed << std::setprecision(6) << result.jaccard << ","
             << std::fixed << std::setprecision(3) << result.theta_degrees << ","
             << std::fixed << std::setprecision(3) << result.phi_degrees << ","
             << std::fixed << std::setprecision(3) << result.alpha_degrees << ","
             << result.true_positive << ","
             << result.true_negative << ","
             << result.false_positive << ","
             << result.false_negative << ","
             << std::fixed << std::setprecision(6) << result.f1_score << ","
             << std::fixed << std::setprecision(6) << result.precision << ","
             << std::fixed << std::setprecision(6) << result.recall << "\n";
    }
    
    file.close();
    std::cout << "Benchmark results saved to " << filename << "\n";
}

// Get list of folders in current directory
std::vector<std::string> get_folders_in_current_directory() {
    std::vector<std::string> folders;
    try {
        for (const auto& entry : fs::directory_iterator(fs::current_path())) {
            if (entry.is_directory()) {
                folders.push_back(entry.path().filename().string());
            }
        }
        std::sort(folders.begin(), folders.end());
    } catch (const std::exception& e) {
        std::cerr << "Error reading directories: " << e.what() << "\n";
    }
    return folders;
}

// Extract folder name from path (e.g., "data3d_1-5" from "data3d_1-5/file.csv")
std::string extract_folder_name(const std::string& folder_path) {
    fs::path p(folder_path);
    return p.filename().string();
}

void rotate_triangle(Triangle& triangle, double theta_degrees, double phi_degrees, double alpha_degrees) {
    if (std::abs(alpha_degrees) < 1e-6) {
        std::cout << "Warning: alpha_degrees is near zero (" << alpha_degrees << "), no rotation applied.\n";
        return;
    }

    // Преобразование углов в радианы
    double theta = theta_degrees * M_PI / 180.0;
    double phi = phi_degrees * M_PI / 180.0;
    double alpha = alpha_degrees * M_PI / 180.0;

    // Вектор оси вращения
    double ux = std::sin(theta) * std::cos(phi);
    double uy = std::sin(theta) * std::sin(phi);
    double uz = std::cos(theta);

    // Матрица поворота (формула Родрига)
    double cos_a = std::cos(alpha);
    double sin_a = std::sin(alpha);
    double one_minus_cos_a = 1.0 - cos_a;
    double R[3][3] = {
        {cos_a + ux * ux * one_minus_cos_a, ux * uy * one_minus_cos_a - uz * sin_a, ux * uz * one_minus_cos_a + uy * sin_a},
        {uy * ux * one_minus_cos_a + uz * sin_a, cos_a + uy * uy * one_minus_cos_a, uy * uz * one_minus_cos_a - ux * sin_a},
        {uz * ux * one_minus_cos_a - uy * sin_a, uz * uy * one_minus_cos_a + ux * sin_a, cos_a + uz * uz * one_minus_cos_a}
    };

    // Вращение точки
    auto rotate_point = [&](double& x, double& y, double& z) {
        double x_new = R[0][0] * x + R[0][1] * y + R[0][2] * z;
        double y_new = R[1][0] * x + R[1][1] * y + R[1][2] * z;
        double z_new = R[2][0] * x + R[2][1] * y + R[2][2] * z;
        x = x_new;
        y = y_new;
        z = z_new;
        };

    // Вращение вершин
    rotate_point(triangle.v1_x, triangle.v1_y, triangle.v1_z);
    rotate_point(triangle.v2_x, triangle.v2_y, triangle.v2_z);
    rotate_point(triangle.v3_x, triangle.v3_y, triangle.v3_z);

    // Вращение нормали
    double normal_x = triangle.normal_x;
    double normal_y = triangle.normal_y;
    double normal_z = triangle.normal_z;
    rotate_point(normal_x, normal_y, normal_z);

    // Нормализация нормали
    double norm = std::sqrt(normal_x * normal_x + normal_y * normal_y + normal_z * normal_z);
    if (norm > 1e-6) {
        triangle.normal_x = normal_x / norm;
        triangle.normal_y = normal_y / norm;
        triangle.normal_z = normal_z / norm;
    }
    else {
        std::cout << "Warning: Normal vector norm is near zero (" << norm << "), skipping normalization.\n";
    }
}

void run_interactive_mode(const std::string& folder_path, int num_files, bool use_threshold, double threshold, bool use_augmentation, double theta_degrees, double phi_degrees, double alpha_degrees) {
    SatelliteDataset dataset(folder_path);
    const auto& file_names = dataset.get_file_names();
    int files_to_process = std::min(num_files, static_cast<int>(file_names.size()));

    // Выбор метода
    std::cout << "\nSelect method:\n1. Pairwise\n2. Grid Revised\n3. Sweep Line\n4. Shadow Map\n5. Uniform Grid\n6. Ray Casting\n0. Exit\nEnter choice (0-6): ";
    int method_choice;
    std::cin >> method_choice;
    if (method_choice == 0) {
        return;
    }
    if (method_choice < 1 || method_choice > 6) {
        std::cout << "Invalid method choice. Choose 0-6.\n";
        return;
    }

    std::string method_name;
    std::function<std::vector<int>(const std::vector<Triangle>&, const std::vector<double>&)> method_func;
    switch (method_choice) {
    case 1:
        method_name = "pairwise";
        method_func = calculate_labels_pairwise;
        break;
    case 2:
        method_name = "zsorted_proj";
        method_func = [use_threshold, threshold](const std::vector<Triangle>& t, const std::vector<double>& s) {
            return calculate_labels_zsorted_proj(t, s);
            };
        break;
    case 3:
        method_name = "sweep_line";
        method_func = calculate_labels_sweep_line;
        break;
    case 4:
        method_name = "shadow_map";
        method_func = calculate_labels_shadow_map;
        break;
    case 5:
        method_name = "uniform_grid";
        method_func = calculate_labels_uniform_grid;
        break;
    case 6:
        method_name = "ray_casting";
        method_func = calculate_labels_ray_casting;
        break;
    }

    std::vector<double> sun_vector = { 0.0, 0.0, 1.0 };

    // Генерация случайных углов для аугментации
    double theta = theta_degrees;
    double phi = phi_degrees;
    double alpha = alpha_degrees;
    if (use_augmentation && theta_degrees == -1.0) {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        std::uniform_real_distribution<> theta_dist(0.0, 180.0);
        std::uniform_real_distribution<> phi_dist(0.0, 360.0);
        std::uniform_real_distribution<> alpha_dist(1.0, 360.0);
        theta = theta_dist(gen);
        phi = phi_dist(gen);
        alpha = alpha_dist(gen);
        std::cout << "Random rotation: theta = " << theta << " degrees, phi = " << phi << " degrees, alpha = " << alpha << " degrees\n";
    }

    // Обработка каждого файла
    for (int file_idx = 0; file_idx < files_to_process; ++file_idx) {
        // Загрузка данных
        dataset.load_single_file(file_names[file_idx]);
        auto triangles_original = dataset.get_triangles();
        std::cout << "\nProcessing file: " << file_names[file_idx] << "\n";
        std::cout << "Loaded " << triangles_original.size() << " triangles\n";

        // Compute reference labels using pairwise method
        std::vector<int> reference_labels;
        if (!use_augmentation) {
            // No augmentation: compute reference on original data
            reference_labels = calculate_labels_pairwise(triangles_original, sun_vector);
        } else {
            // With augmentation: compute reference on rotated data
            auto triangles_ref = triangles_original;
            for (auto& triangle : triangles_ref) {
                rotate_triangle(triangle, theta, phi, alpha);
            }
            reference_labels = calculate_labels_pairwise(triangles_ref, sun_vector);
        }

        // Create working copy of triangles
        auto triangles = triangles_original;

        // Аугментация
        if (use_augmentation) {
            for (auto& triangle : triangles) {
                rotate_triangle(triangle, theta, phi, alpha);
            }
        }

        // Вычисление меток
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<int> labels = method_func(triangles, sun_vector);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        std::cout << "Using " << method_name << " method\n";
        std::cout << "Time: " << duration << " ms\n";

        // Calculate and display metrics
        if (method_name != "pairwise" && labels.size() == reference_labels.size()) {
            double accuracy = calculate_accuracy(labels, reference_labels);
            double dice = calculate_dice(labels, reference_labels);
            double jaccard = calculate_jaccard(labels, reference_labels);
            ConfusionMetrics cm = calculate_confusion_metrics(labels, reference_labels);
            
            std::cout << "\n--- Accuracy Metrics ---\n";
            std::cout << "Accuracy: " << std::fixed << std::setprecision(2) << accuracy << "%\n";
            std::cout << "Dice coefficient: " << std::fixed << std::setprecision(2) << dice << "%\n";
            std::cout << "Jaccard index: " << std::fixed << std::setprecision(2) << jaccard << "%\n";
            std::cout << "F1 Score: " << std::fixed << std::setprecision(4) << cm.f1_score << "\n";
            std::cout << "Precision: " << std::fixed << std::setprecision(4) << cm.precision << "\n";
            std::cout << "Recall: " << std::fixed << std::setprecision(4) << cm.recall << "\n";
            std::cout << "\nConfusion Matrix:\n";
            std::cout << "True Positives: " << cm.true_positive << "\n";
            std::cout << "True Negatives: " << cm.true_negative << "\n";
            std::cout << "False Positives: " << cm.false_positive << "\n";
            std::cout << "False Negatives: " << cm.false_negative << "\n";
        } else if (method_name == "pairwise") {
            std::cout << "\n--- Pairwise Method (Reference) ---\n";
            std::cout << "Accuracy: 100.00% (reference method)\n";
        }

        // Вывод результатов
        std::cout << "\nLabels (first 10): ";
        for (size_t i = 0; i < std::min<size_t>(10, labels.size()); ++i) {
            std::cout << labels[i] << " ";
        }
        std::cout << "\nNumber of shadowed triangles: " << std::count(labels.begin(), labels.end(), 0) << "\n";

        // Подготовка треугольников для визуализации
        for (size_t i = 0; i < triangles.size(); ++i) {
            triangles[i].label = static_cast<double>(labels[i]);
        }

        // Сохранение результатов
        std::string output_filename = method_name + "_" + std::to_string(file_idx + 1) + ".csv";
        save_results(triangles, labels, reference_labels, output_filename);
        std::cout << "Results saved to results_3d/" << output_filename << "\n";

        // Визуализация
        visualize_triangles("results_3d/" + output_filename, 1, false);
    }
}

void run_benchmark_mode(const std::string& folder_path, int num_files, bool use_threshold, double threshold, bool use_augmentation, double theta_degrees, double phi_degrees, double alpha_degrees, bool random_per_run = false) {
    SatelliteDataset dataset(folder_path);
    const auto& file_names = dataset.get_file_names();
    int files_to_process = std::min(num_files, static_cast<int>(file_names.size()));

    const int num_runs = 10;
    std::vector<std::pair<std::string, std::function<std::vector<int>(const std::vector<Triangle>&, const std::vector<double>&)>>> methods = {
        {"pairwise", calculate_labels_pairwise},
        {"zsorted_proj", calculate_labels_zsorted_proj},
        {"sweep_line", calculate_labels_sweep_line},
        {"shadow_map", calculate_labels_shadow_map},
        {"uniform_grid", calculate_labels_uniform_grid},
        {"ray_casting", calculate_labels_ray_casting}
    };
    std::vector<double> sun_vector = { 0.0, 0.0, 1.0 };

    // Vector to store all benchmark results
    std::vector<BenchmarkResult> all_results;

    // Determine augmentation mode: 0 = none, 1 = manual, 2 = random
    int aug_mode = 0;
    if (use_augmentation) {
        aug_mode = random_per_run ? 2 : 1;
    }

    // Setup random number generator with time-based seed for random mode
    std::mt19937 gen;
    std::uniform_real_distribution<> theta_dist(0.0, 180.0);
    std::uniform_real_distribution<> phi_dist(0.0, 360.0);
    std::uniform_real_distribution<> alpha_dist(1.0, 360.0);
    
    if (random_per_run) {
        // Use time-based seed for true randomness
        auto now = std::chrono::high_resolution_clock::now();
        auto seed = std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();
        gen.seed(static_cast<unsigned int>(seed));
    }

    // Generate single random rotation for non-per-run mode
    double theta = theta_degrees;
    double phi = phi_degrees;
    double alpha = alpha_degrees;
    if (use_augmentation && !random_per_run && theta_degrees == -1.0) {
        theta = theta_dist(gen);
        phi = phi_dist(gen);
        alpha = alpha_dist(gen);
        std::cout << "Random rotation: theta = " << theta << " degrees, phi = " << phi << " degrees, alpha = " << alpha << " degrees\n";
    }

    for (int file_idx = 0; file_idx < files_to_process; ++file_idx) {
        dataset.load_single_file(file_names[file_idx]);
        auto triangles_original = dataset.get_triangles();
        std::cout << "\nProcessing file: " << file_names[file_idx] << "\n";
        std::cout << "Loaded " << triangles_original.size() << " triangles\n";

        // OPTIMIZED: Compute reference labels once per file for non-random augmentation
        std::vector<int> reference_labels_fixed;
        if (!random_per_run) {
            if (!use_augmentation) {
                // No augmentation: compute reference on original data once
                reference_labels_fixed = calculate_labels_pairwise(triangles_original, sun_vector);
                std::cout << "Reference labels computed on original data.\n";
            } else {
                // Fixed augmentation: apply rotation once and compute reference
                auto triangles_rotated = triangles_original;
                for (auto& triangle : triangles_rotated) {
                    rotate_triangle(triangle, theta, phi, alpha);
                }
                reference_labels_fixed = calculate_labels_pairwise(triangles_rotated, sun_vector);
                std::cout << "Reference labels computed on rotated data (fixed rotation).\n";

                triangles_original = triangles_rotated;
            }
        }

        std::cout << "\nRunning benchmark for " << num_runs << " runs per method:\n";
        
        for (const auto& [method_name, method_func] : methods) {
            double total_time = 0.0;
            double total_accuracy = 0.0;
            double total_dice = 0.0;
            double total_jaccard = 0.0;
            
            for (int run = 0; run < num_runs; ++run) {

                if (method_name == "zsorted_proj")
                   break;

                if (method_name == "pairwise" && run==2)
                    break;

                // Generate new random angles for each run if in random_per_run mode
                double run_theta = theta;
                double run_phi = phi;
                double run_alpha = alpha;
                auto triangles = triangles_original;

                
                // OPTIMIZED: For random augmentation, generate new angles and compute reference per run
                std::vector<int> reference_labels;
                if (random_per_run && use_augmentation) {
                    // Generate new angles for this specific run
                    run_theta = theta_dist(gen);
                    run_phi = phi_dist(gen);
                    run_alpha = alpha_dist(gen);
                    
                    // Apply rotation and compute reference labels for this run
                    auto triangles_ref = triangles_original;
                    for (auto& triangle : triangles_ref) {
                        rotate_triangle(triangle, run_theta, run_phi, run_alpha);
                    }
                    reference_labels = calculate_labels_pairwise(triangles_ref, sun_vector);

                    triangles = triangles_ref;
                } else {
                    // Use pre-computed reference labels
                    reference_labels = reference_labels_fixed;
                }
              

                // Измерение времени
                auto start = std::chrono::high_resolution_clock::now();
                std::vector<int> labels = method_func(triangles, sun_vector);
                auto end = std::chrono::high_resolution_clock::now();
                double time_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;

                // Расчет метрик
                ConfusionMetrics cm;
                double accuracy = 0.0;
                double dice = 0.0;
                double jaccard = 0.0;
                
                if (method_name != "pairwise") {
                    if (labels.size() == reference_labels.size()) {
                        cm = calculate_confusion_metrics(labels, reference_labels);
                        accuracy = calculate_accuracy(labels, reference_labels);
                        dice = calculate_dice(labels, reference_labels);
                        jaccard = calculate_jaccard(labels, reference_labels);
                    }
                } else {
                    // Pairwise всегда 100% точности
                    accuracy = 100.0;
                    dice = 100.0;
                    jaccard = 100.0;
                    for (size_t j = 0; j < labels.size(); ++j) {
                        if (labels[j] == 1) cm.true_positive++;
                        else cm.true_negative++;
                    }
                    cm.f1_score = 1.0;
                    cm.precision = 1.0;
                    cm.recall = 1.0;
                }
                
                // Store result for this run
                BenchmarkResult result;
                result.method_name = method_name;
                result.file_name = file_names[file_idx];
                result.run_number = run + 1;
                result.num_triangles = static_cast<int>(triangles.size());
                result.time_ms = time_ms;
                result.accuracy = accuracy;
                result.dice = dice;
                result.jaccard = jaccard;
                result.theta_degrees = use_augmentation ? run_theta : 0.0;
                result.phi_degrees = use_augmentation ? run_phi : 0.0;
                result.alpha_degrees = use_augmentation ? run_alpha : 0.0;
                result.true_positive = cm.true_positive;
                result.true_negative = cm.true_negative;
                result.false_positive = cm.false_positive;
                result.false_negative = cm.false_negative;
                result.f1_score = cm.f1_score;
                result.precision = cm.precision;
                result.recall = cm.recall;
                all_results.push_back(result);
                
                total_time += time_ms;
                total_accuracy += accuracy;
                total_dice += dice;
                total_jaccard += jaccard;
            }
            
            int nr = num_runs;
            if (method_name == "pairwise")  nr = 2;

            // Print summary for this method
            std::cout << std::left << std::setw(15) << method_name
                << ": avg time = " << std::fixed << std::setprecision(3) << total_time / nr
                << " ms, avg accuracy = " << std::fixed << std::setprecision(2) << total_accuracy / nr << "%"
                << ", avg dice = " << std::fixed << std::setprecision(2) << total_dice / nr << "%"
                << ", avg jaccard = " << std::fixed << std::setprecision(2) << total_jaccard / nr << "%\n";
        }
    }
    
    // Generate output filename with folder name and augmentation mode
    std::string folder_name = extract_folder_name(folder_path);
    std::string output_filename = "benchmark_results_" + folder_name + "_aug" + std::to_string(aug_mode) + ".csv";
    
    // Save all results to CSV
    save_benchmark_results(all_results, output_filename);
}

void run_comprehensive_test_mode() {
    std::string test_folder = "test_overall_data";
    
    // Check if test folder exists
    if (!fs::exists(test_folder) || !fs::is_directory(test_folder)) {
        std::cout << "Error: Folder '" << test_folder << "' does not exist in current directory.\n";
        std::cout << "Please create it and add data folders (e.g., data3d_1-5, data3d_6-10, etc.)\n";
        return;
    }
    
    // Get all subfolders in test_overall_data
    std::vector<std::string> data_folders;
    try {
        for (const auto& entry : fs::directory_iterator(test_folder)) {
            if (entry.is_directory()) {
                data_folders.push_back(entry.path().string());
            }
        }
        std::sort(data_folders.begin(), data_folders.end());
    } catch (const std::exception& e) {
        std::cerr << "Error reading test folder: " << e.what() << "\n";
        return;
    }
    
    if (data_folders.empty()) {
        std::cout << "No data folders found in '" << test_folder << "'.\n";
        return;
    }
    
    std::cout << "\nFound " << data_folders.size() << " data folders:\n";
    for (size_t i = 0; i < data_folders.size(); ++i) {
        std::cout << "  " << (i + 1) << ". " << extract_folder_name(data_folders[i]) << "\n";
    }
    
    // Submenu for augmentation mode selection
    std::cout << "\nComprehensive Test Mode - Select augmentation:\n";
    std::cout << "1. No augmentation (aug0)\n";
    std::cout << "2. Fixed rotation (aug1) - manual angles\n";
    std::cout << "3. Random rotation (aug2) - changes each run\n";
    std::cout << "4. All three modes\n";
    std::cout << "0. Cancel\n";
    std::cout << "Enter choice (0-4): ";
    
    int aug_choice;
    std::cin >> aug_choice;
    
    if (aug_choice == 0) {
        std::cout << "Cancelled.\n";
        return;
    }
    
    if (aug_choice < 1 || aug_choice > 4) {
        std::cout << "Invalid choice.\n";
        return;
    }
    
    // Test parameters
    int num_files = 100000; // Process all files
    bool use_threshold = false;
    double threshold = 1e-6;
    
    // Manual rotation angles for mode 1
    double manual_theta = 60.0;
    double manual_phi = 0.0;
    double manual_alpha = 135.0;
    
    std::cout << "\nStarting comprehensive test...\n\n";
    
    // Test each folder with selected mode(s)
    for (const auto& folder_path : data_folders) {
        std::string folder_name = extract_folder_name(folder_path);
        std::cout << "\n========================================\n";
        std::cout << "Testing folder: " << folder_name << "\n";
        std::cout << "========================================\n";
        
        // Mode 0: No augmentation
        if (aug_choice == 1 || aug_choice == 4) {
            std::cout << "\n--- Mode 0: No augmentation ---\n";
            run_benchmark_mode(folder_path, num_files, use_threshold, threshold, false, 0.0, 0.0, 90.0, false);
        }
        
        // Mode 1: Manual/Fixed rotation
        if (aug_choice == 2 || aug_choice == 4) {
            std::cout << "\n--- Mode 1: Fixed rotation (theta=" << manual_theta << ", phi=" << manual_phi << ", alpha=" << manual_alpha << ") ---\n";
            run_benchmark_mode(folder_path, num_files, use_threshold, threshold, true, manual_theta, manual_phi, manual_alpha, false);
        }
        
        // Mode 2: Random rotation (per-run)
        if (aug_choice == 3 || aug_choice == 4) {
            std::cout << "\n--- Mode 2: Random rotation (changes each run) ---\n";
            run_benchmark_mode(folder_path, num_files, use_threshold, threshold, true, -1.0, -1.0, -1.0, true);
        }
    }
    
    std::cout << "\n========================================\n";
    std::cout << "Comprehensive test completed!\n";
    std::cout << "========================================\n";
}

int main() {
    try {
        // Параметры
        std::string folder_path = "test";
        int num_files = 1;
        bool use_threshold = false;
        double threshold = 1e-6;
        bool use_augmentation = false;
        double theta_degrees = 0.0;
        double phi_degrees = 0.0;
        double alpha_degrees = 90.0;

        while (true) {
            std::cout << "\nMain Menu:\n"
                << "1. Interactive mode\n"
                << "2. Benchmark mode\n"
                << "3. Comprehensive test mode\n"
                << "4. Set threshold (current: " << (use_threshold ? "ON" : "OFF") << ", threshold: " << threshold << ")\n"
                << "5. Set number of files (current: " << num_files << ")\n"
                << "6. Select data folder (current: " << folder_path << ")\n"
                << "7. Set data augmentation (current: " << (use_augmentation ? "ON" : "OFF") << ", theta: " << theta_degrees << " deg, phi: " << phi_degrees << " deg, alpha: " << alpha_degrees << " deg)\n"
                << "0. Exit\n"
                << "Enter choice (0-7): ";
            int choice;
            std::cin >> choice;

            if (choice == 0) {
                std::cout << "Exiting program.\n";
                break;
            }

            switch (choice) {
            case 1:
                run_interactive_mode(folder_path, num_files, use_threshold, threshold, use_augmentation, theta_degrees, phi_degrees, alpha_degrees);
                break;
            case 2:
                run_benchmark_mode(folder_path, num_files, use_threshold, threshold, use_augmentation, theta_degrees, phi_degrees, alpha_degrees, false);
                break;
            case 3:
                run_comprehensive_test_mode();
                break;
            case 4:
                std::cout << "\nThreshold Settings:\n"
                    << "1. Enable/Disable threshold\n"
                    << "2. Set threshold value\n"
                    << "0. Back\n"
                    << "Enter choice (0-2): ";
                int threshold_choice;
                std::cin >> threshold_choice;
                if (threshold_choice == 1) {
                    std::cout << "Enable threshold? (0 = No, 1 = Yes): ";
                    int enable;
                    std::cin >> enable;
                    use_threshold = (enable == 1);
                }
                else if (threshold_choice == 2) {
                    std::cout << "Enter threshold value (e.g., 1e-6): ";
                    std::cin >> threshold;
                }
                else if (threshold_choice != 0) {
                    std::cout << "Invalid choice. Choose 0-2.\n";
                }
                break;
            case 5:
                std::cout << "Enter number of files: ";
                std::cin >> num_files;
                if (num_files <= 0) {
                    std::cout << "Number of files must be positive.\n";
                    num_files = 1;
                }
                break;
            case 6: {
                std::cout << "\nAvailable folders in current directory:\n";
                auto folders = get_folders_in_current_directory();
                if (folders.empty()) {
                    std::cout << "No folders found.\n";
                    break;
                }
                for (size_t i = 0; i < folders.size(); ++i) {
                    std::cout << "  " << (i + 1) << ". " << folders[i] << "\n";
                }
                std::cout << "Enter folder number (or 0 to cancel): ";
                int folder_choice;
                std::cin >> folder_choice;
                if (folder_choice > 0 && folder_choice <= static_cast<int>(folders.size())) {
                    folder_path = folders[folder_choice - 1];
                    std::cout << "Selected folder: " << folder_path << "\n";
                } else if (folder_choice != 0) {
                    std::cout << "Invalid folder selection.\n";
                }
                break;
            }
            case 7:
                std::cout << "\nAugmentation Settings:\n"
                    << "1. Enable/Disable augmentation\n"
                    << "2. Random rotation\n"
                    << "3. Manual rotation\n"
                    << "0. Back\n"
                    << "Enter choice (0-3): ";
                int aug_choice;
                std::cin >> aug_choice;
                if (aug_choice == 1) {
                    std::cout << "Enable augmentation? (0 = No, 1 = Yes): ";
                    int enable;
                    std::cin >> enable;
                    use_augmentation = (enable == 1);
                }
                else if (aug_choice == 2) {
                    use_augmentation = true;
                    theta_degrees = -1.0;
                    phi_degrees = -1.0;
                    alpha_degrees = -1.0;
                    std::cout << "Random rotation enabled.\n";
                }
                else if (aug_choice == 3) {
                    use_augmentation = true;
                    std::cout << "Enter theta (latitude, 0-180 degrees): ";
                    std::cin >> theta_degrees;
                    std::cout << "Enter phi (longitude, 0-360 degrees): ";
                    std::cin >> phi_degrees;
                    do {
                        std::cout << "Enter alpha (rotation angle, degrees, non-zero): ";
                        std::cin >> alpha_degrees;
                        if (std::abs(alpha_degrees) < 1e-6) {
                            std::cout << "Error: alpha must be non-zero. Try again.\n";
                        }
                    } while (std::abs(alpha_degrees) < 1e-6);
                }
                else if (aug_choice != 0) {
                    std::cout << "Invalid choice. Choose 0-3.\n";
                }
                break;
            default:
                std::cout << "Invalid choice. Choose 0-7.\n";
            }
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
