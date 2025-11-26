#pragma once

#include <vector>
#include <string>
#include "SatelliteDataset.h"

std::vector<int> calculate_labels_pairwise(const std::vector<Triangle>& triangles, const std::vector<double>& sun_vector);

std::vector<int> calculate_labels_zsorted_proj(const std::vector<Triangle>& triangles, const std::vector<double>& sun_vector6);

std::vector<int> calculate_labels_sweep_line(const std::vector<Triangle>& triangles, const std::vector<double>& sun_vector);

std::vector<int> calculate_labels_shadow_map(const std::vector<Triangle>& triangles, const std::vector<double>& sun_vector);

std::vector<int> calculate_labels_uniform_grid(const std::vector<Triangle>& triangles, const std::vector<double>& sun_vector);

std::vector<int> calculate_labels_ray_casting(const std::vector<Triangle>& triangles, const std::vector<double>& sun_vector);



void visualize_triangles(const std::string& csv_file, int every, bool show_normals);

void save_results(const std::vector<Triangle>& triangles, const std::vector<int>& labels, const std::vector<int>& reference_labels, const std::string& filename);
