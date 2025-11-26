#include "ShadowAlgorithms.h"
#include <fstream>
#include <stdexcept>
#include <cmath>
#include <pybind11/embed.h>
#include <pybind11/stl.h>
#include <filesystem>
namespace py = pybind11;


namespace fs = std::filesystem;

static struct PyInitializer {
    PyInitializer() {
        PyConfig config;
        PyConfig_InitPythonConfig(&config);

        // Helper function to check status and throw exceptions
        auto check_status = [](PyStatus status, const std::string& msg) {
            if (PyStatus_IsError(status)) {
                throw std::runtime_error(msg + ": " + status.err_msg);
            }
            };

        // Configure Python paths
        const std::wstring python_path = L"C:\\Users\\zapevalin\\AppData\\Local\\Programs\\Python\\Python312";
        const std::wstring python_exe = python_path + L"\\python.exe";
        const std::wstring lib_path = python_path + L"\\Lib";
        const std::wstring site_packages_path = python_path + L"\\Lib\\site-packages";
        const std::wstring dlls_path = python_path + L"\\DLLs";
        const std::wstring project_path = fs::current_path().wstring();

        // Configure program_name, executable and home
        check_status(PyConfig_SetString(&config, &config.program_name, python_exe.c_str()),
            "Failed to set program_name");
        check_status(PyConfig_SetString(&config, &config.executable, python_exe.c_str()),
            "Failed to set executable");
        check_status(PyConfig_SetString(&config, &config.home, python_path.c_str()),
            "Failed to set PYTHONHOME");

        // Add paths to module_search_paths
        config.module_search_paths_set = 1;
        const std::wstring* paths[] = {
            &lib_path,
            &site_packages_path,
            &dlls_path,
            &project_path
        };

        for (const auto* path : paths) {
            check_status(PyWideStringList_Append(&config.module_search_paths, path->c_str()),
                "Failed to append path");
        }

        // Initialize Python interpreter
        check_status(Py_InitializeFromConfig(&config), "Failed to initialize Python");
        PyConfig_Clear(&config);
    }

    ~PyInitializer() { Py_Finalize(); }
} py_init;




#include "ShadowAlgorithms.h"
#include <fstream>
#include <stdexcept>
#include <cmath>
#include <set>
#include <pybind11/embed.h>
#include <pybind11/stl.h>
#include <filesystem>
#include <numeric> // for std::iota
#include <iostream> // for std::cerr

#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>

#include <vector>
#include <array>
#include <algorithm>
#include <limits>
#include <stdexcept>
namespace py = pybind11;


/*TRY*/
#include <vector>
#include <array>
#include <algorithm>
#include <limits>
#include <stdexcept>





// Structure for AABB (Axis-Aligned Bounding Box)
struct AABB {
    double min_x, min_y, min_z;
    double max_x, max_y, max_z;
};

// BVH Node structure
struct BVHNode {
    AABB bounds;
    size_t triangle_idx; // Triangle index if leaf node
    size_t left, right;  // Indices of child nodes if not leaf
    bool is_leaf;
};



bool ray_triangle_intersection(
    const double* origin, const double* dir,
    const double* v0, const double* v1, const double* v2,
    double& t) {
    const double EPSILON = 1e-6;
    double edge1[3], edge2[3], h[3], s[3], q[3];
    double a, f, u, v;

    // Compute edge vectors
    for (int i = 0; i < 3; ++i) {
        edge1[i] = v1[i] - v0[i];
        edge2[i] = v2[i] - v0[i];
    }

    // h = cross(dir, edge2)
    h[0] = dir[1] * edge2[2] - dir[2] * edge2[1];
    h[1] = dir[2] * edge2[0] - dir[0] * edge2[2];
    h[2] = dir[0] * edge2[1] - dir[1] * edge2[0];

    // a = dot(edge1, h)
    a = edge1[0] * h[0] + edge1[1] * h[1] + edge1[2] * h[2];

    if (std::abs(a) < EPSILON) {
        return false; // Ray parallel to triangle
    }

    f = 1.0 / a;
    for (int i = 0; i < 3; ++i) {
        s[i] = origin[i] - v0[i];
    }

    // u = f * dot(s, h)
    u = f * (s[0] * h[0] + s[1] * h[1] + s[2] * h[2]);
    if (u < 0.0 || u > 1.0) {
        return false;
    }

    // q = cross(s, edge1)
    q[0] = s[1] * edge1[2] - s[2] * edge1[1];
    q[1] = s[2] * edge1[0] - s[0] * edge1[2];
    q[2] = s[0] * edge1[1] - s[1] * edge1[0];

    // v = f * dot(dir, q)
    v = f * (dir[0] * q[0] + dir[1] * q[1] + dir[2] * q[2]);
    if (v < 0.0 || u + v > 1.0) {
        return false;
    }

    // t = f * dot(edge2, q)
    t = f * (edge2[0] * q[0] + edge2[1] * q[1] + edge2[2] * q[2]);
    if (t > EPSILON) {
        return true;
    }

    return false;
}



std::vector<int> calculate_labels_pairwise(const std::vector<Triangle>& triangles, const std::vector<double>& sun_vector) {
    if (triangles.empty()) {
        return std::vector<int>();
    }
    if (sun_vector.size() != 3) {
        throw std::runtime_error("Sun vector must have 3 components");
    }

    size_t N = triangles.size();
    const double epsilon = 1e-6; // Tolerance for floating-point comparisons

    std::vector<int> labels(N, 1);

    std::vector<double> sca(N);
    for (size_t i = 0; i < N; ++i) {
        sca[i] = triangles[i].normal_x * sun_vector[0] +
            triangles[i].normal_y * sun_vector[1] +
            triangles[i].normal_z * sun_vector[2];
        if (sca[i] <= epsilon) {
            labels[i] = 0;
        }
    }

    for (size_t i = 0; i < N; ++i) {
        if (labels[i] == 0) continue;

        double center[3] = {
            (triangles[i].v1_x + triangles[i].v2_x + triangles[i].v3_x) / 3.0,
            (triangles[i].v1_y + triangles[i].v2_y + triangles[i].v3_y) / 3.0,
            (triangles[i].v1_z + triangles[i].v2_z + triangles[i].v3_z) / 3.0
        };

        for (size_t j = 0; j < N; ++j) {
            if (i == j || sca[j] <= epsilon) continue;

            double v0[3] = { triangles[j].v1_x, triangles[j].v1_y, triangles[j].v1_z };
            double v1[3] = { triangles[j].v2_x, triangles[j].v2_y, triangles[j].v2_z };
            double v2[3] = { triangles[j].v3_x, triangles[j].v3_y, triangles[j].v3_z };

            double center_j[3] = {
                (triangles[j].v1_x + triangles[j].v2_x + triangles[j].v3_x) / 3.0,
                (triangles[j].v1_y + triangles[j].v2_y + triangles[j].v3_y) / 3.0,
                (triangles[j].v1_z + triangles[j].v2_z + triangles[j].v3_z) / 3.0
            };
            double proj_i = center[0] * sun_vector[0] + center[1] * sun_vector[1] + center[2] * sun_vector[2];
            double proj_j = center_j[0] * sun_vector[0] + center_j[1] * sun_vector[1] + center_j[2] * sun_vector[2];
            if (proj_j < proj_i - epsilon) continue;

            double t;
            if (ray_triangle_intersection(center, sun_vector.data(), v0, v1, v2, t) && t >= -epsilon) {
                labels[i] = 0;
                break;
            }
        }
    }

    return labels;
}


std::vector<int> calculate_labels_ray_casting(const std::vector<Triangle>& triangles, const std::vector<double>& sun_vector) {
    if (triangles.empty()) {
        return std::vector<int>();
    }
    if (sun_vector.size() != 3 || sun_vector[0] != 0.0 || sun_vector[1] != 0.0 || sun_vector[2] != 1.0) {
        throw std::runtime_error("Sun vector must be {0, 0, 1}");
    }

    size_t N = triangles.size();
    const double epsilon = 1e-6;

    // Initialize labels and centers
    std::vector<int> labels(N, 1);
    std::vector<double> sca(N);
    std::vector<std::array<double, 3>> centers(N);
    for (size_t i = 0; i < N; ++i) {
        centers[i] = {
            (triangles[i].v1_x + triangles[i].v2_x + triangles[i].v3_x) / 3.0,
            (triangles[i].v1_y + triangles[i].v2_y + triangles[i].v3_y) / 3.0,
            (triangles[i].v1_z + triangles[i].v2_z + triangles[i].v3_z) / 3.0
        };
        sca[i] = triangles[i].normal_z;
        if (sca[i] <= epsilon) {
            labels[i] = 0;
        }
    }

    // Build BVH
    std::vector<AABB> triangle_bounds(N);
    for (size_t i = 0; i < N; ++i) {
        triangle_bounds[i] = {
            std::min({triangles[i].v1_x, triangles[i].v2_x, triangles[i].v3_x}),
            std::min({triangles[i].v1_y, triangles[i].v2_y, triangles[i].v3_y}),
            std::min({triangles[i].v1_z, triangles[i].v2_z, triangles[i].v3_z}),
            std::max({triangles[i].v1_x, triangles[i].v2_x, triangles[i].v3_x}),
            std::max({triangles[i].v1_y, triangles[i].v2_y, triangles[i].v3_y}),
            std::max({triangles[i].v1_z, triangles[i].v2_z, triangles[i].v3_z})
        };
    }

    std::vector<BVHNode> nodes;
    std::vector<size_t> indices(N);
    std::iota(indices.begin(), indices.end(), 0);

    // Recursive BVH construction
    auto build_bvh = [&](auto& self, size_t start, size_t end, size_t depth) -> size_t {
        BVHNode node;
        if (end - start == 1) {
            node.is_leaf = true;
            node.triangle_idx = indices[start];
            node.bounds = triangle_bounds[node.triangle_idx];
            nodes.push_back(node);
            return nodes.size() - 1;
        }

        // Choose axis for split (by largest extent)
        double x_range = 0, y_range = 0, z_range = 0;
        double x_min = std::numeric_limits<double>::infinity(), x_max = -x_min;
        double y_min = x_min, y_max = -x_min;
        double z_min = x_min, z_max = -x_min;
        for (size_t i = start; i < end; ++i) {
            x_min = std::min(x_min, triangle_bounds[indices[i]].min_x);
            x_max = std::max(x_max, triangle_bounds[indices[i]].max_x);
            y_min = std::min(y_min, triangle_bounds[indices[i]].min_y);
            y_max = std::max(y_max, triangle_bounds[indices[i]].max_y);
            z_min = std::min(z_min, triangle_bounds[indices[i]].min_z);
            z_max = std::max(z_max, triangle_bounds[indices[i]].max_z);
        }
        x_range = x_max - x_min;
        y_range = y_max - y_min;
        z_range = z_max - z_min;

        int axis = 0;
        if (y_range > x_range && y_range > z_range) axis = 1;
        else if (z_range > x_range) axis = 2;

        // Sort along chosen axis
        std::sort(indices.begin() + start, indices.begin() + end, [&](size_t a, size_t b) {
            return centers[a][axis] < centers[b][axis];
        });

        size_t mid = (start + end) / 2;
        node.is_leaf = false;
        node.left = self(self, start, mid, depth + 1);
        node.right = self(self, mid, end, depth + 1);

        // Merge bounds
        node.bounds = {
            std::min(nodes[node.left].bounds.min_x, nodes[node.right].bounds.min_x),
            std::min(nodes[node.left].bounds.min_y, nodes[node.right].bounds.min_y),
            std::min(nodes[node.left].bounds.min_z, nodes[node.right].bounds.min_z),
            std::max(nodes[node.left].bounds.max_x, nodes[node.right].bounds.max_x),
            std::max(nodes[node.left].bounds.max_y, nodes[node.right].bounds.max_y),
            std::max(nodes[node.left].bounds.max_z, nodes[node.right].bounds.max_z)
        };
        nodes.push_back(node);
        return nodes.size() - 1;
    };

    size_t root = build_bvh(build_bvh, 0, N, 0);

    // Check ray-AABB intersection
    auto ray_aabb_intersection = [](const double origin[3], const double dir[3], const AABB& box, double& t_min, double& t_max) {
        double t1 = (box.min_x - origin[0]) / (dir[0] + 1e-10);
        double t2 = (box.max_x - origin[0]) / (dir[0] + 1e-10);
        t_min = std::min(t1, t2);
        t_max = std::max(t1, t2);
        t1 = (box.min_y - origin[1]) / (dir[1] + 1e-10);
        t2 = (box.max_y - origin[1]) / (dir[1] + 1e-10);
        t_min = std::max(t_min, std::min(t1, t2));
        t_max = std::min(t_max, std::max(t1, t2));
        t1 = (box.min_z - origin[2]) / (dir[2] + 1e-10);
        t2 = (box.max_z - origin[2]) / (dir[2] + 1e-10);
        t_min = std::max(t_min, std::min(t1, t2));
        t_max = std::min(t_max, std::max(t1, t2));
        return t_max >= t_min && t_max >= -1e-6;
    };

    // Check for shadowing
    for (size_t i = 0; i < N; ++i) {
        if (labels[i] == 0) continue;
        double origin[3] = { centers[i][0], centers[i][1], centers[i][2] };
        double dir[3] = { sun_vector[0], sun_vector[1], sun_vector[2] };
        std::vector<size_t> stack = { root };
        bool shadowed = false;

        while (!stack.empty() && !shadowed) {
            size_t node_idx = stack.back();
            stack.pop_back();
            const BVHNode& node = nodes[node_idx];

            double t_min, t_max;
            if (!ray_aabb_intersection(origin, dir, node.bounds, t_min, t_max)) continue;

            if (node.is_leaf) {
                size_t j = node.triangle_idx;
                if (i == j || sca[j] <= epsilon || centers[j][2] <= centers[i][2] + epsilon) continue;
                double v0[3] = { triangles[j].v1_x, triangles[j].v1_y, triangles[j].v1_z };
                double v1[3] = { triangles[j].v2_x, triangles[j].v2_y, triangles[j].v2_z };
                double v2[3] = { triangles[j].v3_x, triangles[j].v3_y, triangles[j].v3_z };
                double t;
                if (ray_triangle_intersection(origin, dir, v0, v1, v2, t)) {
                    labels[i] = 0;
                    shadowed = true;
                }
            } else {
                stack.push_back(node.left);
                stack.push_back(node.right);
            }
        }
    }

    return labels;
}

std::vector<int> calculate_labels_uniform_grid(const std::vector<Triangle>& triangles, const std::vector<double>& sun_vector) {
    if (triangles.empty()) {
        return std::vector<int>();
    }
    if (sun_vector.size() != 3 || sun_vector[0] != 0.0 || sun_vector[1] != 0.0 || sun_vector[2] != 1.0) {
        throw std::runtime_error("Sun vector must be {0, 0, 1}");
    }

    size_t N = triangles.size();
    const double epsilon = 1e-6;

    // Initialize labels and centers
    std::vector<int> labels(N, 1);
    std::vector<double> sca(N);
    std::vector<std::array<double, 3>> centers(N);
    for (size_t i = 0; i < N; ++i) {
        centers[i] = {
            (triangles[i].v1_x + triangles[i].v2_x + triangles[i].v3_x) / 3.0,
            (triangles[i].v1_y + triangles[i].v2_y + triangles[i].v3_y) / 3.0,
            (triangles[i].v1_z + triangles[i].v2_z + triangles[i].v3_z) / 3.0
        };
        sca[i] = triangles[i].normal_z;
        if (sca[i] <= epsilon) {
            labels[i] = 0;
        }
    }

    // Determine mesh bounds
    double x_min = std::numeric_limits<double>::infinity(), x_max = -x_min;
    double y_min = x_min, y_max = -x_min;
    double z_min = x_min, z_max = -x_min;
    for (size_t i = 0; i < N; ++i) {
        x_min = std::min(x_min, std::min({ triangles[i].v1_x, triangles[i].v2_x, triangles[i].v3_x }));
        x_max = std::max(x_max, std::max({ triangles[i].v1_x, triangles[i].v2_x, triangles[i].v3_x }));
        y_min = std::min(y_min, std::min({ triangles[i].v1_y, triangles[i].v2_y, triangles[i].v3_y }));
        y_max = std::max(y_max, std::max({ triangles[i].v1_y, triangles[i].v2_y, triangles[i].v3_y }));
        z_min = std::min(z_min, std::min({ triangles[i].v1_z, triangles[i].v2_z, triangles[i].v3_z }));
        z_max = std::max(z_max, std::max({ triangles[i].v1_z, triangles[i].v2_z, triangles[i].v3_z }));
    }

    // Create 3D grid
    int num_bins = std::max(10, static_cast<int>(std::cbrt(N) * 2.0)); // Adaptive bin count for 3D
    double x_bin_width = (x_max - x_min) / num_bins;
    double y_bin_width = (y_max - y_min) / num_bins;
    double z_bin_width = (z_max - z_min) / num_bins;

    // Grid to store potential occluders
    std::vector<std::vector<std::vector<std::vector<size_t>>>> grid(
        num_bins, std::vector<std::vector<std::vector<size_t>>>(
            num_bins, std::vector<std::vector<size_t>>(
                num_bins)));

    // Insert triangles into grid
    for (size_t i = 0; i < N; ++i) {
        double x_min_tri = std::min({ triangles[i].v1_x, triangles[i].v2_x, triangles[i].v3_x });
        double x_max_tri = std::max({ triangles[i].v1_x, triangles[i].v2_x, triangles[i].v3_x });
        double y_min_tri = std::min({ triangles[i].v1_y, triangles[i].v2_y, triangles[i].v3_y });
        double y_max_tri = std::max({ triangles[i].v1_y, triangles[i].v2_y, triangles[i].v3_y });
        double z_min_tri = std::min({ triangles[i].v1_z, triangles[i].v2_z, triangles[i].v3_z });
        double z_max_tri = std::max({ triangles[i].v1_z, triangles[i].v2_z, triangles[i].v3_z });

        int x_start = std::max(0, static_cast<int>((x_min_tri - x_min) / x_bin_width));
        int x_end = std::min(num_bins, static_cast<int>((x_max_tri - x_min) / x_bin_width) + 1);
        int y_start = std::max(0, static_cast<int>((y_min_tri - y_min) / y_bin_width));
        int y_end = std::min(num_bins, static_cast<int>((y_max_tri - y_min) / y_bin_width) + 1);
        int z_start = std::max(0, static_cast<int>((z_min_tri - z_min) / z_bin_width));
        int z_end = std::min(num_bins, static_cast<int>((z_max_tri - z_min) / z_bin_width) + 1);

        for (int x = x_start; x < x_end; ++x) {
            for (int y = y_start; y < y_end; ++y) {
                for (int z = z_start; z < z_end; ++z) {
                    grid[x][y][z].push_back(i);
                }
            }
        }
    }

    // Check point in triangle (2D, for projection on XY plane)
    auto point_in_triangle = [epsilon](double px, double py, double x1, double y1, double x2, double y2, double x3, double y3, double& u, double& v) -> bool {
        double v0x = x2 - x1, v0y = y2 - y1;
        double v1x = x3 - x1, v1y = y3 - y1;
        double v2x = px - x1, v2y = py - y1;

        double dot00 = v0x * v0x + v0y * v0y;
        double dot01 = v0x * v1x + v0y * v1y;
        double dot02 = v0x * v2x + v0y * v2y;
        double dot11 = v1x * v1x + v1y * v1y;
        double dot12 = v1x * v2x + v1y * v2y;

        double denom = dot00 * dot11 - dot01 * dot01;
        if (std::abs(denom) < epsilon) {
            u = v = 0.0;
            return false;
        }
        double invDenom = 1.0 / denom;
        u = (dot11 * dot02 - dot01 * dot12) * invDenom;
        v = (dot00 * dot12 - dot01 * dot02) * invDenom;

        return (u >= -epsilon) && (v >= -epsilon) && (u + v <= 1.0 + epsilon);
        };

    // Check for shadowing
    const double threshold_normal = 0.1;
    for (size_t i = 0; i < N; ++i) {
        if (labels[i] == 0) continue;
        double x_center = centers[i][0];
        double y_center = centers[i][1];
        double z_center = centers[i][2];
        int x_idx = std::max(0, std::min(num_bins - 1, static_cast<int>((x_center - x_min) / x_bin_width)));
        int y_idx = std::max(0, std::min(num_bins - 1, static_cast<int>((y_center - y_min) / y_bin_width)));
        int z_idx = std::max(0, std::min(num_bins - 1, static_cast<int>((z_center - z_min) / z_bin_width)));

        bool shadowed = false;
        for (int dx = -1; dx <= 1 && !shadowed; ++dx) {
            for (int dy = -1; dy <= 1 && !shadowed; ++dy) {
                for (int dz = 0; dz <= num_bins && !shadowed; ++dz) { // Check only cells above in Z
                    int nx = x_idx + dx, ny = y_idx + dy, nz = z_idx + dz;
                    if (nx < 0 || nx >= num_bins || ny < 0 || ny >= num_bins || nz >= num_bins) continue;
                    for (size_t j : grid[nx][ny][nz]) {
                        if (i == j || sca[j] <= epsilon) continue;
                        if (centers[j][2] <= z_center + epsilon) continue; // Check only triangles above
                        double x1 = triangles[j].v1_x, y1 = triangles[j].v1_y, z1 = triangles[j].v1_z;
                        double x2 = triangles[j].v2_x, y2 = triangles[j].v2_y, z2 = triangles[j].v2_z;
                        double x3 = triangles[j].v3_x, y3 = triangles[j].v3_y, z3 = triangles[j].v3_z;
                        double u, v;
                        if (point_in_triangle(x_center, y_center, x1, y1, x2, y2, x3, y3, u, v)) {
                            double interpolated_z = z1 * (1 - u - v) + z2 * u + z3 * v;
                            if (z_center < interpolated_z - epsilon) {
                                shadowed = true;
                                labels[i] = 0;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    return labels;
}


std::vector<int> calculate_labels_zsorted_proj(const std::vector<Triangle>& triangles, const std::vector<double>& sun_vector) {
    if (triangles.empty()) {
        return std::vector<int>();
    }
    if (sun_vector.size() != 3 || sun_vector[0] != 0.0 || sun_vector[1] != 0.0 || sun_vector[2] != 1.0) {
        throw std::runtime_error("Sun vector must be {0, 0, 1}");
    }

    size_t N = triangles.size();
    const double epsilon = 1e-6;

    // Precompute data (for optimization)
    std::vector<std::array<double, 3>> centers(N);
    std::vector<double> normal_z(N);
    std::vector<int> labels(N, 1);
    std::vector<double> sca(N);

    for (size_t i = 0; i < N; ++i) {
        centers[i] = {
            (triangles[i].v1_x + triangles[i].v2_x + triangles[i].v3_x) / 3.0,
            (triangles[i].v1_y + triangles[i].v2_y + triangles[i].v3_y) / 3.0,
            (triangles[i].v1_z + triangles[i].v2_z + triangles[i].v3_z) / 3.0
        };
        normal_z[i] = triangles[i].normal_z;
        sca[i] = normal_z[i];
        if (sca[i] <= epsilon) {
            labels[i] = 0;
        }
    }

    // Check if point is in triangle with epsilon tolerance and return barycentric coordinates for Z interpolation
    auto point_in_triangle_bary = [epsilon](double px, double py, double x1, double y1, double x2, double y2, double x3, double y3, double& u, double& v) -> bool {
        double v0x = x2 - x1, v0y = y2 - y1;
        double v1x = x3 - x1, v1y = y3 - y1;
        double v2x = px - x1, v2y = py - y1;

        double dot00 = v0x * v0x + v0y * v0y;
        double dot01 = v0x * v1x + v0y * v1y;
        double dot02 = v0x * v2x + v0y * v2y;
        double dot11 = v1x * v1x + v1y * v1y;
        double dot12 = v1x * v2x + v1y * v2y;

        double denom = dot00 * dot11 - dot01 * dot01;
        if (std::abs(denom) < epsilon) {
            u = v = 0.0;
            return false;
        }
        double invDenom = 1.0 / denom;
        u = (dot11 * dot02 - dot01 * dot12) * invDenom;
        v = (dot00 * dot12 - dot01 * dot02) * invDenom;

        return (u >= -epsilon) && (v >= -epsilon) && (u + v <= 1.0 + epsilon);
        };

    // Sort by Z (for optimization)
    std::vector<size_t> sort_indices(N);
    std::iota(sort_indices.begin(), sort_indices.end(), 0);
    std::sort(sort_indices.begin(), sort_indices.end(), [&](size_t a, size_t b) {
        if (std::abs(centers[a][2] - centers[b][2]) > epsilon) {
            return centers[a][2] < centers[b][2];
        }
        return (centers[a][0] < centers[b][0] || (centers[a][0] == centers[b][0] && centers[a][1] < centers[b][1]));
        });

    // Reorder arrays (for optimization)
    std::vector<std::array<double, 3>> sorted_centers(N);
    std::vector<double> sorted_sca(N);
    std::vector<int> sorted_labels(N);
    for (size_t i = 0; i < N; ++i) {
        sorted_centers[i] = centers[sort_indices[i]];
        sorted_sca[i] = sca[sort_indices[i]];
        sorted_labels[i] = labels[sort_indices[i]];
    }

    // Check shadowing with reverse Z-order and barycentric interpolation
    const double threshold_normal = 0.1;
    double comp_threshold =  epsilon;
    for (size_t i = 0; i < N; ++i) {
        if (sorted_labels[i] == 0 || std::abs(sorted_sca[i]) < threshold_normal) continue;
        double x1 = triangles[sort_indices[i]].v1_x, y1 = triangles[sort_indices[i]].v1_y, z1 = triangles[sort_indices[i]].v1_z;
        double x2 = triangles[sort_indices[i]].v2_x, y2 = triangles[sort_indices[i]].v2_y, z2 = triangles[sort_indices[i]].v2_z;
        double x3 = triangles[sort_indices[i]].v3_x, y3 = triangles[sort_indices[i]].v3_y, z3 = triangles[sort_indices[i]].v3_z;

        for (size_t j = 0; j < i; ++j) {
            if (sorted_labels[j] == 0 || sorted_sca[j] <= epsilon || std::abs(sorted_sca[j]) < threshold_normal) continue;
            double x_center_j = sorted_centers[j][0];
            double y_center_j = sorted_centers[j][1];
            double z_center_j = sorted_centers[j][2];

            double u, v;
            if (point_in_triangle_bary(x_center_j, y_center_j, x1, y1, x2, y2, x3, y3, u, v)) {
                // Interpolate Z at centroid location on triangle i
                double interpolated_z = z1 * (1 - u - v) + z2 * u + z3 * v;
                if (z_center_j < interpolated_z - comp_threshold) {
                    sorted_labels[j] = 0;
                }
            }
        }
    }

    // Additional check for nearly vertical triangles (if still lit after first pass)
    for (size_t i = 0; i < N; ++i) {
        if (sorted_labels[i] == 0 || std::abs(sorted_sca[i]) >= threshold_normal) continue;
        double z_center_i = sorted_centers[i][2];
        bool shadowed = false;
        for (size_t j = i + 1; j < N; ++j) {  // Check forward
            if (std::abs(sorted_sca[j]) < threshold_normal) continue;
            double x1_j = triangles[sort_indices[j]].v1_x, y1_j = triangles[sort_indices[j]].v1_y, z1_j = triangles[sort_indices[j]].v1_z;
            double x2_j = triangles[sort_indices[j]].v2_x, y2_j = triangles[sort_indices[j]].v2_y, z2_j = triangles[sort_indices[j]].v2_z;
            double x3_j = triangles[sort_indices[j]].v3_x, y3_j = triangles[sort_indices[j]].v3_y, z3_j = triangles[sort_indices[j]].v3_z;

            double u, v;
            if (point_in_triangle_bary(sorted_centers[i][0], sorted_centers[i][1], x1_j, y1_j, x2_j, y2_j, x3_j, y3_j, u, v)) {
                double interpolated_z_j = z1_j * (1 - u - v) + z2_j * u + z3_j * v;
                if (z_center_i < interpolated_z_j - comp_threshold) {
                    shadowed = true;
                    break;
                }
            }
        }
        if (shadowed) sorted_labels[i] = 0;
    }

    // Restore original order
    std::vector<int> result_labels(N);
    for (size_t i = 0; i < N; ++i) {
        result_labels[sort_indices[i]] = sorted_labels[i];
    }

    return result_labels;
}


std::vector<int> calculate_labels_sweep_line(const std::vector<Triangle>& triangles, const std::vector<double>& sun_vector) {
    if (triangles.empty()) {
        return std::vector<int>();
    }
    if (sun_vector.size() != 3 || sun_vector[0] != 0.0 || sun_vector[1] != 0.0 || sun_vector[2] != 1.0) {
        throw std::runtime_error("Sun vector must be {0, 0, 1}");
    }

    size_t N = triangles.size();
    const double epsilon = 1e-6;

    // Precompute data
    std::vector<std::array<double, 3>> centers(N);
    std::vector<double> normal_z(N);
    std::vector<int> labels(N, 1);
    std::vector<double> sca(N);

    for (size_t i = 0; i < N; ++i) {
        centers[i] = {
            (triangles[i].v1_x + triangles[i].v2_x + triangles[i].v3_x) / 3.0,
            (triangles[i].v1_y + triangles[i].v2_y + triangles[i].v3_y) / 3.0,
            (triangles[i].v1_z + triangles[i].v2_z + triangles[i].v3_z) / 3.0
        };
        normal_z[i] = triangles[i].normal_z;
        sca[i] = normal_z[i];
        if (sca[i] <= epsilon) {
            labels[i] = 0;
        }
    }

    // Check point in triangle with barycentric coordinates
    auto point_in_triangle = [epsilon](double px, double py, double x1, double y1, double x2, double y2, double x3, double y3, double& u, double& v) -> bool {
        double v0x = x2 - x1, v0y = y2 - y1;
        double v1x = x3 - x1, v1y = y3 - y1;
        double v2x = px - x1, v2y = py - y1;

        double dot00 = v0x * v0x + v0y * v0y;
        double dot01 = v0x * v1x + v0y * v1y;
        double dot02 = v0x * v2x + v0y * v2y;
        double dot11 = v1x * v1x + v1y * v1y;
        double dot12 = v1x * v2x + v1y * v2y;

        double denom = dot00 * dot11 - dot01 * dot01;
        if (std::abs(denom) < epsilon) {
            u = v = 0.0;
            return false;
        }
        double invDenom = 1.0 / denom;
        u = (dot11 * dot02 - dot01 * dot12) * invDenom;
        v = (dot00 * dot12 - dot01 * dot02) * invDenom;

        return (u >= -epsilon) && (v >= -epsilon) && (u + v <= 1.0 + epsilon);
        };

    // Create events
    std::vector<std::tuple<double, double, size_t, int>> events;
    events.reserve(2 * N);
    for (size_t i = 0; i < N; ++i) {
        events.emplace_back(centers[i][2], centers[i][1], i, 0); // check event
        if (sca[i] > epsilon && std::abs(sca[i]) >= 0.1) {
            events.emplace_back(centers[i][2], centers[i][1], i, 1); // start event
        }
    }

    // Sort by Z, Y, type
    std::sort(events.begin(), events.end(), [](const auto& a, const auto& b) {
        if (std::abs(std::get<0>(a) - std::get<0>(b)) > 1e-6) {
            return std::get<0>(a) > std::get<0>(b);
        }
        if (std::abs(std::get<1>(a) - std::get<1>(b)) > 1e-6) {
            return std::get<1>(a) > std::get<1>(b);
        }
        return std::get<3>(a) < std::get<3>(b);
        });

    // Active triangles
    struct ActiveTriangle {
        double x_center, y_center;
        size_t idx;
        bool operator<(const ActiveTriangle& other) const {
            if (std::abs(x_center - other.x_center) > 1e-6) return x_center < other.x_center;
            return y_center < other.y_center;
        }
    };
    std::set<ActiveTriangle> active_triangles;

    // Process events
    double comp_threshold = epsilon;
    const double threshold_normal = 0.1;
    for (const auto& [z, y, idx, event_type] : events) {
        // Remove obsolete triangles
        auto it = active_triangles.begin();
        while (it != active_triangles.end() && centers[it->idx][2] < z - epsilon) {
            it = active_triangles.erase(it);
        }

        if (event_type == 1) { // start event
            active_triangles.insert({ centers[idx][0], centers[idx][1], idx });
        }
        else { // check event
            if (labels[idx] == 0) continue;
            double x_center = centers[idx][0];
            double y_center = centers[idx][1];
            double z_center = centers[idx][2];
            bool shadowed = false;
            for (const auto& tri : active_triangles) {
                if (std::abs(sca[tri.idx]) < threshold_normal) continue; // Skip vertical triangles
                double x1 = triangles[tri.idx].v1_x, y1 = triangles[tri.idx].v1_y, z1 = triangles[tri.idx].v1_z;
                double x2 = triangles[tri.idx].v2_x, y2 = triangles[tri.idx].v2_y, z2 = triangles[tri.idx].v2_z;
                double x3 = triangles[tri.idx].v3_x, y3 = triangles[tri.idx].v3_y, z3 = triangles[tri.idx].v3_z;
                double u, v;
                if (point_in_triangle(x_center, y_center, x1, y1, x2, y2, x3, y3, u, v)) {
                    double interpolated_z = z1 * (1 - u - v) + z2 * u + z3 * v;
                    if (z_center < interpolated_z - comp_threshold) {
                        shadowed = true;
                        break;
                    }
                }
            }
            if (shadowed) labels[idx] = 0;
            // Additional check for vertical triangles at check event
            if (!shadowed && std::abs(sca[idx]) < threshold_normal) {
                for (const auto& tri : active_triangles) {
                    if (std::abs(sca[tri.idx]) < threshold_normal) continue;
                    double x1 = triangles[tri.idx].v1_x, y1 = triangles[tri.idx].v1_y, z1 = triangles[tri.idx].v1_z;
                    double x2 = triangles[tri.idx].v2_x, y2 = triangles[tri.idx].v2_y, z2 = triangles[tri.idx].v2_z;
                    double x3 = triangles[tri.idx].v3_x, y3 = triangles[tri.idx].v3_y, z3 = triangles[tri.idx].v3_z;
                    double u, v;
                    if (point_in_triangle(x_center, y_center, x1, y1, x2, y2, x3, y3, u, v)) {
                        double interpolated_z = z1 * (1 - u - v) + z2 * u + z3 * v;
                        if (z_center < interpolated_z - comp_threshold) {
                            shadowed = true;
                            break;
                        }
                    }
                }
                if (shadowed) labels[idx] = 0;
            }
        }
    }

    return labels;
}


std::vector<int> calculate_labels_shadow_map(const std::vector<Triangle>& triangles, const std::vector<double>& sun_vector) {
    if (triangles.empty()) {
        return std::vector<int>();
    }
    if (sun_vector.size() != 3 || sun_vector[0] != 0.0 || sun_vector[1] != 0.0 || sun_vector[2] != 1.0) {
        throw std::runtime_error("Sun vector must be {0, 0, 1}");
    }

    size_t N = triangles.size();
    const double epsilon = 1e-6; // Tolerance for floating-point comparisons

    // Precompute data
    struct TriangleData {
        double center[3];
        double normal_z;
    };
    std::vector<TriangleData> tri_data(N);
    std::vector<int> labels(N, 1);
    std::vector<double> sca(N);

    for (size_t i = 0; i < N; ++i) {
        tri_data[i].center[0] = (triangles[i].v1_x + triangles[i].v2_x + triangles[i].v3_x) / 3.0;
        tri_data[i].center[1] = (triangles[i].v1_y + triangles[i].v2_y + triangles[i].v3_y) / 3.0;
        tri_data[i].center[2] = (triangles[i].v1_z + triangles[i].v2_z + triangles[i].v3_z) / 3.0;
        tri_data[i].normal_z = triangles[i].normal_z;
        sca[i] = tri_data[i].normal_z;
        if (sca[i] <= epsilon) {
            labels[i] = 0;
        }
    }

    // Determine X and Y bounds for shadow map
    double x_min = std::numeric_limits<double>::infinity(), x_max = -std::numeric_limits<double>::infinity();
    double y_min = std::numeric_limits<double>::infinity(), y_max = -std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < N; ++i) {
        x_min = std::min(x_min, std::min({ triangles[i].v1_x, triangles[i].v2_x, triangles[i].v3_x }));
        x_max = std::max(x_max, std::max({ triangles[i].v1_x, triangles[i].v2_x, triangles[i].v3_x }));
        y_min = std::min(y_min, std::min({ triangles[i].v1_y, triangles[i].v2_y, triangles[i].v3_y }));
        y_max = std::max(y_max, std::max({ triangles[i].v1_y, triangles[i].v2_y, triangles[i].v3_y }));
    }

    // Create 2D grid
    int num_bins_x = std::max(100, static_cast<int>(std::sqrt(N) * 2.20));
    int num_bins_y = std::max(100, static_cast<int>(std::sqrt(N) * 2.20));
    std::vector<double> x_edges(num_bins_x + 1);
    std::vector<double> y_edges(num_bins_y + 1);
    double x_bin_width = (x_max - x_min) / num_bins_x;
    double y_bin_width = (y_max - y_min) / num_bins_y;
    for (int i = 0; i <= num_bins_x; ++i) {
        x_edges[i] = x_min + i * x_bin_width;
    }
    for (int i = 0; i <= num_bins_y; ++i) {
        y_edges[i] = y_min + i * y_bin_width;
    }

    // Initialize shadow map
    std::vector<std::vector<double>> shadow_map(num_bins_x, std::vector<double>(num_bins_y, -std::numeric_limits<double>::infinity()));
    std::vector<std::vector<std::vector<size_t>>> triangle_indices(num_bins_x, std::vector<std::vector<size_t>>(num_bins_y));

    // Optimized function to check point in triangle
    auto point_in_triangle = [epsilon](double px, double py, double x1, double y1, double x2, double y2, double x3, double y3) -> bool {
        auto cross = [](double x1, double y1, double x2, double y2, double px, double py) {
            return (px - x2) * (y1 - y2) - (x1 - x2) * (py - y2);
            };
        double d1 = cross(px, py, x1, y1, x2, y2);
        double d2 = cross(px, py, x2, y2, x3, y3);
        double d3 = cross(px, py, x3, y3, x1, y1);
        bool has_neg = (d1 < -epsilon) || (d2 < -epsilon) || (d3 < -epsilon);
        bool has_pos = (d1 > epsilon) || (d2 > epsilon) || (d3 > epsilon);
        return !(has_neg && has_pos);
        };

    // Populate shadow map, skipping triangles with |normal_z| < threshold_normal
    const double threshold_normal = 0.1;
    for (size_t i = 0; i < N; ++i) {
        if (sca[i] <= epsilon || std::abs(sca[i]) < threshold_normal) continue; // Skip back-facing and vertical triangles
        double x1 = triangles[i].v1_x, y1 = triangles[i].v1_y;
        double x2 = triangles[i].v2_x, y2 = triangles[i].v2_y;
        double x3 = triangles[i].v3_x, y3 = triangles[i].v3_y;
        double z_center = tri_data[i].center[2];

        // Conservative rasterization
        double x_min_tri = std::min({ x1, x2, x3 });
        double x_max_tri = std::max({ x1, x2, x3 });
        double y_min_tri = std::min({ y1, y2, y3 });
        double y_max_tri = std::max({ y1, y2, y3 });
        int x_start = static_cast<int>(std::floor((x_min_tri - x_min) / x_bin_width));
        int x_end = static_cast<int>(std::ceil((x_max_tri - x_min) / x_bin_width));
        int y_start = static_cast<int>(std::floor((y_min_tri - y_min) / y_bin_width));
        int y_end = static_cast<int>(std::ceil((y_max_tri - y_min) / y_bin_width));
        x_start = std::max(0, x_start);
        x_end = std::min(num_bins_x, x_end);
        y_start = std::max(0, y_start);
        y_end = std::min(num_bins_y, y_end);

        // Rasterize bins
        for (int x_idx = x_start; x_idx < x_end; ++x_idx) {
            for (int y_idx = y_start; y_idx < y_end; ++y_idx) {
                double x_bin = x_edges[x_idx] + x_bin_width * 0.5;
                double y_bin = y_edges[y_idx] + y_bin_width * 0.5;
                if (point_in_triangle(x_bin, y_bin, x1, y1, x2, y2, x3, y3)) {
                    shadow_map[x_idx][y_idx] = std::max(shadow_map[x_idx][y_idx], z_center);
                    triangle_indices[x_idx][y_idx].push_back(i);
                }
            }
        }
    }

    // Check triangles against shadow map with threshold
    double comp_threshold =  epsilon;
    for (size_t i = 0; i < N; ++i) {
        if (labels[i] == 0) continue;
        double x_center = tri_data[i].center[0];
        double y_center = tri_data[i].center[1];
        double z_center = tri_data[i].center[2];
        int x_idx = static_cast<int>(std::floor((x_center - x_min) / x_bin_width));
        int y_idx = static_cast<int>(std::floor((y_center - y_min) / y_bin_width));
        x_idx = std::max(0, std::min(x_idx, num_bins_x - 1));
        y_idx = std::max(0, std::min(y_idx, num_bins_y - 1));

        // For triangles with |normal_z| < threshold_normal use detailed check
        if (std::abs(sca[i]) < threshold_normal) {
            double max_z = -std::numeric_limits<double>::infinity();
            for (size_t idx : triangle_indices[x_idx][y_idx]) {
                if (idx != i) {
                    max_z = std::max(max_z, tri_data[idx].center[2]);
                }
            }
            if (z_center < max_z - comp_threshold) {
                labels[i] = 0;
            }
        }
        else {
            // Standard check
            if (z_center < shadow_map[x_idx][y_idx] - comp_threshold) {
                labels[i] = 0;
            }
        }
    }

    return labels;
}



void visualize_triangles(const std::string& csv_file, int every, bool show_normals) {
    try {
        py::module_ visualize = py::module_::import("visualize3d");
        visualize.attr("visualize_triangles")(csv_file, every, show_normals);
    }
    catch (const py::error_already_set& e) {
        std::cerr << "Python error: " << e.what() << "\n";
        throw;
    }
}

void save_results(const std::vector<Triangle>& triangles, const std::vector<int>& labels, const std::vector<int>& reference_labels, const std::string& filename) {
    std::filesystem::create_directory("results_3d");
    std::string full_path = "results_3d/" + filename;

    std::ofstream file(full_path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open " + full_path + " for writing");
    }

    file << "Triangle ID,Component Type,V1_X,V1_Y,V1_Z,V2_X,V2_Y,V2_Z,V3_X,V3_Y,V3_Z,Normal_X,Normal_Y,Normal_Z,Label,Correct\n";

    for (size_t i = 0; i < triangles.size(); ++i) {
        const auto& t = triangles[i];
        int correct = (labels[i] == reference_labels[i]) ? 1 : 0;
        file << i << "," << t.component_type << ","
            << t.v1_x << "," << t.v1_y << "," << t.v1_z << ","
            << t.v2_x << "," << t.v2_y << "," << t.v2_z << ","
            << t.v3_x << "," << t.v3_y << "," << t.v3_z << ","
            << t.normal_x << "," << t.normal_y << "," << t.normal_z << ","
            << labels[i] << "," << correct << "\n";
    }

    file.close();
}
