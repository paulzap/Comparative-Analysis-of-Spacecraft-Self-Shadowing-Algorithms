#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <numbers>
#include <random>
#include <iostream>

const int BODY_min = 5;
const int BODY_max = 10;
const int SOLAR_min = 5;
const int SOLAR_max = 15;
const int ANTENNA_min = 1;
const int ANTENNA_max = 5;
const int OFFSET = 3;

struct Point {
    double x, y, z;
    Point operator+(const Point& other) const { return { x + other.x, y + other.y, z + other.z }; }
    Point operator-(const Point& other) const { return { x - other.x, y - other.y, z - other.z }; }
    Point operator*(double scalar) const { return { x * scalar, y * scalar, z * scalar }; }
};

struct Triangle {
    Point v1, v2, v3;
};

class Shape {
protected:
    std::string name;
    Point center;
    std::vector<Point> vertices;
    std::vector<double> parameters;
    double polygon_size;

public:
    const void setPolygonSize(double pol_size) { polygon_size = pol_size; }
    virtual void generateBodyVerticies() {};
    virtual void generateComponentVerticies(const std::vector<Point>& _limits, int _place) = 0;
    virtual void generateRandomParameters() {}
    virtual void generateRandomParametersSolar() {}
    virtual void generateRandomParametersAntenna() {}
    virtual std::vector<Triangle> getTriangles() const = 0;
    virtual void storeParameters(std::ofstream& out) const = 0;
    virtual std::vector<Point> computeTriangleNormals(const Triangle& triangle) const = 0;
    virtual void generateCustomParameters(const std::vector<double>& _params) {}
    virtual const int getSolarClass() const { return -1; }
    virtual ~Shape() = default;

    const std::vector<Point>& getVertices() const { return vertices; }
    const std::vector<double>& getParameters() const { return parameters; }
    const Point& getCenter() const { return center; }

    virtual std::vector<Point> getExtremePoints() const {
        if (vertices.empty()) {
            throw std::runtime_error("Vertices vector is empty");
        }
        Point left = vertices[0], right = vertices[0], bottom = vertices[0], top = vertices[0], front = vertices[0], back = vertices[0];
        for (const auto& vertex : vertices) {
            if (vertex.x < left.x) left = vertex;
            if (vertex.x > right.x) right = vertex;
            if (vertex.y < bottom.y) bottom = vertex;
            if (vertex.y > top.y) top = vertex;
            if (vertex.z < front.z) front = vertex;
            if (vertex.z > back.z) back = vertex;
        }
        return { right, bottom, left, top, front, back };
    }

    double generateRandomDouble(int min, int max) {
        return (double)(rand() % (max - min + 1) + min);
    }

protected:
    virtual std::vector<Triangle> triangulateFaceTriangles(const std::vector<Point>& face, double polygon_size, bool reverse = false) const {
        std::vector<Triangle> triangles;
        if (face.size() != 4) {
            if (reverse) {
                triangles.push_back({ face[0], face[2], face[1] });
            }
            else {
                triangles.push_back({ face[0], face[1], face[2] });
            }
            return triangles;
        }

        double min_x = face[0].x, max_x = face[0].x;
        double min_y = face[0].y, max_y = face[0].y;
        double min_z = face[0].z, max_z = face[0].z;
        for (const auto& p : face) {
            min_x = std::min(min_x, p.x);
            max_x = std::max(max_x, p.x);
            min_y = std::min(min_y, p.y);
            max_y = std::max(max_y, p.y);
            min_z = std::min(min_z, p.z);
            max_z = std::max(max_z, p.z);
        }

        double dx = max_x - min_x;
        double dy = max_y - min_y;
        double dz = max_z - min_z;
        bool xy_plane = dz < dx && dz < dy;
        bool xz_plane = dy < dx && dy < dz;
        bool yz_plane = dx < dy && dx < dz;

        int nx = xy_plane ? static_cast<int>(std::ceil(dx / polygon_size)) : xz_plane ? static_cast<int>(std::ceil(dx / polygon_size)) : static_cast<int>(std::ceil(dy / polygon_size));
        int ny = xy_plane ? static_cast<int>(std::ceil(dy / polygon_size)) : xz_plane ? static_cast<int>(std::ceil(dz / polygon_size)) : static_cast<int>(std::ceil(dz / polygon_size));
        if (nx < 1) nx = 1;
        if (ny < 1) ny = 1;

        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                double x0, x1, y0, y1, z0, z1;
                if (xy_plane) {
                    x0 = min_x + i * polygon_size;
                    x1 = min_x + (i + 1) * polygon_size;
                    y0 = min_y + j * polygon_size;
                    y1 = min_y + (j + 1) * polygon_size;
                    z0 = min_z;
                    if (x1 > max_x) x1 = max_x;
                    if (y1 > max_y) y1 = max_y;
                }
                else if (xz_plane) {
                    x0 = min_x + i * polygon_size;
                    x1 = min_x + (i + 1) * polygon_size;
                    z0 = min_z + j * polygon_size;
                    z1 = min_z + (j + 1) * polygon_size;
                    y0 = min_y;
                    if (x1 > max_x) x1 = max_x;
                    if (z1 > max_z) z1 = max_z;
                }
                else {
                    y0 = min_y + i * polygon_size;
                    y1 = min_y + (i + 1) * polygon_size;
                    z0 = min_z + j * polygon_size;
                    z1 = min_z + (j + 1) * polygon_size;
                    x0 = min_x;
                    if (y1 > max_y) y1 = max_y;
                    if (z1 > max_z) z1 = max_z;
                }

                std::vector<Point> tri1, tri2;
                if (xy_plane) {
                    tri1 = { {x0, y0, z0}, {x1, y0, z0}, {x0, y1, z0} };
                    tri2 = { {x1, y0, z0}, {x1, y1, z0}, {x0, y1, z0} };
                }
                else if (xz_plane) {
                    tri1 = { {x0, y0, z0}, {x1, y0, z0}, {x0, y0, z1} };
                    tri2 = { {x1, y0, z0}, {x1, y0, z1}, {x0, y0, z1} };
                }
                else {
                    tri1 = { {x0, y0, z0}, {x0, y1, z0}, {x0, y0, z1} };
                    tri2 = { {x0, y1, z0}, {x0, y1, z1}, {x0, y0, z1} };
                }
                if (reverse) {
                    triangles.push_back({ tri1[0], tri1[2], tri1[1] });
                    triangles.push_back({ tri2[0], tri2[2], tri2[1] });
                }
                else {
                    triangles.push_back({ tri1[0], tri1[1], tri1[2] });
                    triangles.push_back({ tri2[0], tri2[1], tri2[2] });
                }
            }
        }
        return triangles;
    }
};



class Cube : public Shape {
    double side;
public:
    void generateRandomParameters() override {
        name = "cube";
        parameters.clear();
        side = generateRandomDouble(BODY_min, BODY_max);
        parameters.push_back(side);
    }

    void generateRandomParametersSolar() override {
        name = "cube";
        parameters.clear();
        side = generateRandomDouble(SOLAR_min, SOLAR_max);
        parameters.push_back(side);
    }

    void generateBodyVerticies() override {
        vertices.clear();
        center = { 0, 0, 0 };
        double halfSide = side / 2;
        vertices = {
            {-halfSide, halfSide, halfSide}, {halfSide, halfSide, halfSide}, {halfSide, -halfSide, halfSide}, {-halfSide, -halfSide, halfSide},
            {-halfSide, halfSide, -halfSide}, {halfSide, halfSide, -halfSide}, {halfSide, -halfSide, -halfSide}, {-halfSide, -halfSide, -halfSide}
        };
    }

    void generateCustomParameters(const std::vector<double>& _params) override {
        name = "cube";
        parameters.clear();
        side = _params[0];
        parameters.push_back(side);
    }

    const int getSolarClass() const override { return 1; }

    void generateComponentVerticies(const std::vector<Point>& _limits, int _place) override {
        vertices.clear();
        double halfSide = side / 2;
        switch (_place) {
        case 0: center = { _limits[0].x + halfSide + OFFSET, 0, 0 }; break;
        case 1: center = { 0, _limits[1].y - halfSide - OFFSET, 0 }; break;
        case 2: center = { _limits[2].x - halfSide - OFFSET, 0, 0 }; break;
        case 3: center = { 0, _limits[3].y + halfSide + OFFSET, 0 }; break;
        }
        vertices = {
            {center.x - halfSide, center.y + halfSide, center.z + halfSide},
            {center.x + halfSide, center.y + halfSide, center.z + halfSide},
            {center.x + halfSide, center.y - halfSide, center.z + halfSide},
            {center.x - halfSide, center.y - halfSide, center.z + halfSide},
            {center.x - halfSide, center.y + halfSide, center.z - halfSide},
            {center.x + halfSide, center.y + halfSide, center.z - halfSide},
            {center.x + halfSide, center.y - halfSide, center.z - halfSide},
            {center.x - halfSide, center.y - halfSide, center.z - halfSide}
        };
    }

    std::vector<Triangle> getTriangles() const override {
        std::vector<Triangle> triangles;
        std::vector<std::vector<Point>> faces = {
            {vertices[0], vertices[3], vertices[2], vertices[1]}, // �������� ����� (Z = halfSide)
            {vertices[4], vertices[5], vertices[6], vertices[7]}, // ������ ����� (Z = -halfSide)
            {vertices[0], vertices[1], vertices[5], vertices[4]}, // ������� ����� (Y = halfSide)
            {vertices[2], vertices[3], vertices[7], vertices[6]}, // ������ ����� (Y = -halfSide)
            {vertices[1], vertices[5], vertices[6], vertices[2]}, // ������ ����� (X = halfSide)
            {vertices[0], vertices[4], vertices[7], vertices[3]}  // ����� ����� (X = -halfSide)
        };
        for (const auto& face : faces) {
            auto faceTriangles = triangulateFaceTriangles(face, polygon_size);
            triangles.insert(triangles.end(), faceTriangles.begin(), faceTriangles.end());
        }
        return triangles;
    }

    std::vector<Point> computeTriangleNormals(const Triangle& triangle) const override {
        Point vec1 = { triangle.v2.x - triangle.v1.x, triangle.v2.y - triangle.v1.y, triangle.v2.z - triangle.v1.z };
        Point vec2 = { triangle.v3.x - triangle.v1.x, triangle.v3.y - triangle.v1.y, triangle.v3.z - triangle.v1.z };
        Point normal = {
            vec1.y * vec2.z - vec1.z * vec2.y,
            vec1.z * vec2.x - vec1.x * vec2.z,
            vec1.x * vec2.y - vec1.y * vec2.x
        };
        double length = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
        if (length > 0) {
            normal.x /= length;
            normal.y /= length;
            normal.z /= length;
        }
        Point triangle_center = {
            (triangle.v1.x + triangle.v2.x + triangle.v3.x) / 3.0,
            (triangle.v1.y + triangle.v2.y + triangle.v3.y) / 3.0,
            (triangle.v1.z + triangle.v2.z + triangle.v3.z) / 3.0
        };
        Point vec_to_center = {
            triangle_center.x - center.x,
            triangle_center.y - center.y,
            triangle_center.z - center.z
        };
        double dot_product = normal.x * vec_to_center.x + normal.y * vec_to_center.y + normal.z * vec_to_center.z;
        if (dot_product < 0) {
            normal.x = -normal.x;
            normal.y = -normal.y;
            normal.z = -normal.z;
        }
        return { normal };
    }

    void storeParameters(std::ofstream& out) const override {
        out << "\tSide: " << side << "\n";
    }
};

class Plane : public Shape {
    double width, height;
public:
    void generateRandomParametersSolar() override {
        name = "plane";
        parameters.clear();
        width = generateRandomDouble(SOLAR_min, SOLAR_max);
        height = generateRandomDouble(SOLAR_min, SOLAR_max);
        parameters.push_back(width);
        parameters.push_back(height);
    }

    void generateCustomParameters(const std::vector<double>& _params) override {
        name = "plane";
        parameters.clear();
        width = _params[0];
        height = _params[1];
        parameters.push_back(width);
        parameters.push_back(height);
    }

    const int getSolarClass() const override { return 2; }

    void generateComponentVerticies(const std::vector<Point>& _limits, int _place) override {
        vertices.clear();
        double halfWidth = width / 2, halfHeight = height / 2;
        switch (_place) {
        case 0: center = { _limits[0].x + halfWidth + OFFSET, 0, 0 }; break;
        case 1: center = { 0, _limits[1].y - halfHeight - OFFSET, 0 }; break;
        case 2: center = { _limits[2].x - halfWidth - OFFSET, 0, 0 }; break;
        case 3: center = { 0, _limits[3].y + halfHeight + OFFSET, 0 }; break;
        }
        vertices = {
            {center.x - halfWidth, center.y + halfHeight, center.z},
            {center.x + halfWidth, center.y + halfHeight, center.z},
            {center.x + halfWidth, center.y - halfHeight, center.z},
            {center.x - halfWidth, center.y - halfHeight, center.z}
        };
    }

    std::vector<Triangle> getTriangles() const override {
        std::vector<Triangle> triangles;
        auto faceTriangles = triangulateFaceTriangles(vertices, polygon_size);
        triangles.insert(triangles.end(), faceTriangles.begin(), faceTriangles.end());
        std::vector<Point> reversed = { vertices[0], vertices[3], vertices[2], vertices[1] };
        auto backTriangles = triangulateFaceTriangles(reversed, polygon_size, true);
        triangles.insert(triangles.end(), backTriangles.begin(), backTriangles.end());
        return triangles;
    }

    std::vector<Point> computeTriangleNormals(const Triangle& triangle) const override {
        Point vec1 = { triangle.v2.x - triangle.v1.x, triangle.v2.y - triangle.v1.y, triangle.v2.z - triangle.v1.z };
        Point vec2 = { triangle.v3.x - triangle.v1.x, triangle.v3.y - triangle.v1.y, triangle.v3.z - triangle.v1.z };
        Point normal = {
            vec1.y * vec2.z - vec1.z * vec2.y,
            vec1.z * vec2.x - vec1.x * vec2.z,
            vec1.x * vec2.y - vec1.y * vec2.x
        };
        double length = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
        if (length > 0) {
            normal.x /= length;
            normal.y /= length;
            normal.z /= length;
        }
        return { normal, {-normal.x, -normal.y, -normal.z} };
    }

    void storeParameters(std::ofstream& out) const override {
        out << "\tWidth: " << width << "\n";
        out << "\tHeight: " << height << "\n";
    }
};


class Sphere : public Shape {
    double radius;

public:
    void generateRandomParameters() override {
        name = "sphere";
        parameters.clear();
        radius = generateRandomDouble(BODY_min, BODY_max);
        parameters.push_back(radius);
    }

    void generateRandomParametersAntenna() override {
        name = "sphere";
        parameters.clear();
        radius = generateRandomDouble(ANTENNA_min, ANTENNA_max);
        parameters.push_back(radius);
    }

    void generateBodyVerticies() override {
        vertices.clear();
        center = { 0, 0, 0 };
        int RINGS = std::max(4, static_cast<int>(std::numbers::pi * radius / polygon_size));
        int SECTORS = std::max(4, static_cast<int>(2 * std::numbers::pi * radius / polygon_size));
        RINGS = std::min(RINGS, 100); // ����������� ��� ������� ��������
        SECTORS = std::min(SECTORS, 100);
        for (int i = 0; i <= RINGS; ++i) {
            double theta = i * std::numbers::pi / RINGS;
            for (int j = 0; j <= SECTORS; ++j) {
                double phi = j * 2 * std::numbers::pi / SECTORS;
                vertices.push_back({
                    center.x + radius * std::sin(theta) * std::cos(phi),
                    center.y + radius * std::sin(theta) * std::sin(phi),
                    center.z + radius * std::cos(theta)
                    });
            }
        }
    }

    void generateComponentVerticies(const std::vector<Point>& _limits, int _place) override {
        vertices.clear();
        double halfLength = radius;
        switch (_place) {
        case 0: center = { _limits[0].x + halfLength + OFFSET, 0, 0 }; break;
        case 1: center = { 0, _limits[1].y - halfLength - OFFSET, 0 }; break;
        case 2: center = { _limits[2].x - halfLength - OFFSET, 0, 0 }; break;
        case 3: center = { 0, _limits[3].y + halfLength + OFFSET, 0 }; break;
        }
        int RINGS = std::max(4, static_cast<int>(std::numbers::pi * radius / polygon_size));
        int SECTORS = std::max(4, static_cast<int>(2 * std::numbers::pi * radius / polygon_size));
        RINGS = std::min(RINGS, 100);
        SECTORS = std::min(SECTORS, 100);
        for (int i = 0; i <= RINGS; ++i) {
            double theta = i * std::numbers::pi / RINGS;
            for (int j = 0; j <= SECTORS; ++j) {
                double phi = j * 2 * std::numbers::pi / SECTORS;
                vertices.push_back({
                    center.x + radius * std::sin(theta) * std::cos(phi),
                    center.y + radius * std::sin(theta) * std::sin(phi),
                    center.z + radius * std::cos(theta)
                    });
            }
        }
    }

    std::vector<Point> getExtremePoints() const override {
        return {
            {center.x + radius, center.y, center.z},
            {center.x, center.y - radius, center.z},
            {center.x - radius, center.y, center.z},
            {center.x, center.y + radius, center.z},
            {center.x, center.y, center.z - radius},
            {center.x, center.y, center.z + radius}
        };
    }

    std::vector<Triangle> getTriangles() const override {
        std::vector<Triangle> triangles;
        int RINGS = std::max(4, static_cast<int>(std::numbers::pi * radius / polygon_size));
        int SECTORS = std::max(4, static_cast<int>(2 * std::numbers::pi * radius / polygon_size));
        RINGS = std::min(RINGS, 100);
        SECTORS = std::min(SECTORS, 100);
        for (int i = 0; i < RINGS; ++i) {
            for (int j = 0; j < SECTORS; ++j) {
                int idx1 = i * (SECTORS + 1) + j;
                int idx2 = idx1 + 1;
                int idx3 = (i + 1) * (SECTORS + 1) + j;
                int idx4 = idx3 + 1;
                triangles.push_back({ vertices[idx1], vertices[idx2], vertices[idx3] });
                triangles.push_back({ vertices[idx2], vertices[idx4], vertices[idx3] });
            }
        }
        return triangles;
    }

    std::vector<Point> computeTriangleNormals(const Triangle& triangle) const override {
        Point triangle_center = {
            (triangle.v1.x + triangle.v2.x + triangle.v3.x) / 3.0,
            (triangle.v1.y + triangle.v2.y + triangle.v3.y) / 3.0,
            (triangle.v1.z + triangle.v2.z + triangle.v3.z) / 3.0
        };
        Point normal = {
            triangle_center.x - center.x,
            triangle_center.y - center.y,
            triangle_center.z - center.z
        };
        double length = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
        if (length > 0) {
            normal.x /= length;
            normal.y /= length;
            normal.z /= length;
        }
        return { normal, {-normal.x, -normal.y, -normal.z} };
    }

    void storeParameters(std::ofstream& out) const override {
        out << "\tRadius: " << radius << "\n";
    }
};


class CurvedSurface : public Shape {
    double radius, angle, height;

public:
    void generateRandomParametersAntenna() override {
        name = "curved_surface";
        parameters.clear();
        radius = generateRandomDouble(ANTENNA_min, ANTENNA_max);
        angle = generateRandomDouble(90, 180);
        height = generateRandomDouble(ANTENNA_min, ANTENNA_max);
        parameters.push_back(radius);
        parameters.push_back(angle);
        parameters.push_back(height);
    }

    void generateCustomParameters(const std::vector<double>& _params) override {
        name = "curved_surface";
        parameters.clear();
        radius = _params[0];
        angle = _params[1];
        height = _params[2];
        parameters.push_back(radius);
        parameters.push_back(angle);
        parameters.push_back(height);
    }

    void generateComponentVerticies(const std::vector<Point>& _limits, int _place) override {
        vertices.clear();
        double halfLength = radius;
        switch (_place) {
        case 0: center = { _limits[0].x + halfLength + OFFSET, 0, 0 }; break;
        case 1: center = { 0, _limits[1].y - halfLength - OFFSET, 0 }; break;
        case 2: center = { _limits[2].x - halfLength - OFFSET, 0, 0 }; break;
        case 3: center = { 0, _limits[3].y + halfLength + OFFSET, 0 }; break;
        }
        double angle_rad = angle * std::numbers::pi / 180.0;
        double half_height = height / 2;
        int SEGMENTS = std::max(4, static_cast<int>(radius * angle_rad / polygon_size));
        int HEIGHT_SEGMENTS = std::max(4, static_cast<int>(height / polygon_size));
        SEGMENTS = std::min(SEGMENTS, 100);
        HEIGHT_SEGMENTS = std::min(HEIGHT_SEGMENTS, 100);
        for (int i = 0; i <= SEGMENTS; ++i) {
            double a = -angle_rad / 2 + i * angle_rad / SEGMENTS;
            for (int j = 0; j <= HEIGHT_SEGMENTS; ++j) {
                double z = center.z - half_height + j * height / HEIGHT_SEGMENTS;
                vertices.push_back({
                    center.x + radius * std::cos(a),
                    center.y + radius * std::sin(a),
                    z
                    });
            }
        }
    }

    std::vector<Triangle> getTriangles() const override {
        std::vector<Triangle> triangles;
        int SEGMENTS = std::max(4, static_cast<int>(radius * angle * std::numbers::pi / 180.0 / polygon_size));
        int HEIGHT_SEGMENTS = std::max(4, static_cast<int>(height / polygon_size));
        SEGMENTS = std::min(SEGMENTS, 100);
        HEIGHT_SEGMENTS = std::min(HEIGHT_SEGMENTS, 100);
        for (int i = 0; i < SEGMENTS; ++i) {
            for (int j = 0; j < HEIGHT_SEGMENTS; ++j) {
                int idx1 = i * (HEIGHT_SEGMENTS + 1) + j;
                int idx2 = idx1 + 1;
                int idx3 = (i + 1) * (HEIGHT_SEGMENTS + 1) + j;
                int idx4 = idx3 + 1;
                triangles.push_back({ vertices[idx1], vertices[idx2], vertices[idx3] });
                triangles.push_back({ vertices[idx2], vertices[idx4], vertices[idx3] });
                triangles.push_back({ vertices[idx1], vertices[idx3], vertices[idx2] });
                triangles.push_back({ vertices[idx2], vertices[idx3], vertices[idx4] });
            }
        }
        return triangles;
    }

    std::vector<Point> computeTriangleNormals(const Triangle& triangle) const override {
        Point triangle_center = {
            (triangle.v1.x + triangle.v2.x + triangle.v3.x) / 3.0,
            (triangle.v1.y + triangle.v2.y + triangle.v3.y) / 3.0,
            (triangle.v1.z + triangle.v2.z + triangle.v3.z) / 3.0
        };
        Point normal = {
            triangle_center.x - center.x,
            triangle_center.y - center.y,
            0
        };
        double length = std::sqrt(normal.x * normal.x + normal.y * normal.y);
        
        // If radial vector is too small (degenerate triangle or on axis), use cross product
        if (length < 1e-10) {
            Point vec1 = { triangle.v2.x - triangle.v1.x, triangle.v2.y - triangle.v1.y, triangle.v2.z - triangle.v1.z };
            Point vec2 = { triangle.v3.x - triangle.v1.x, triangle.v3.y - triangle.v1.y, triangle.v3.z - triangle.v1.z };
            normal = {
                vec1.y * vec2.z - vec1.z * vec2.y,
                vec1.z * vec2.x - vec1.x * vec2.z,
                vec1.x * vec2.y - vec1.y * vec2.x
            };
            length = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
        }
        
        if (length > 0) {
            normal.x /= length;
            normal.y /= length;
            normal.z /= length;
        }
        return { normal, {-normal.x, -normal.y, -normal.z} };
    }

    void storeParameters(std::ofstream& out) const override {
        out << "\tRadius: " << radius << "\n";
        out << "\tAngle: " << angle << "\n";
        out << "\tHeight: " << height << "\n";
    }
};

class SphericalCap : public Shape {
    double radius, theta_max;

public:
    void generateRandomParameters() override {
        name = "spherical_cap";
        parameters.clear();
        radius = generateRandomDouble(BODY_min, BODY_max);
        theta_max = generateRandomDouble(30, 80);
        parameters.push_back(radius);
        parameters.push_back(theta_max);
    }

    void generateRandomParametersAntenna() override {
        name = "spherical_cap";
        parameters.clear();
        radius = generateRandomDouble(ANTENNA_min, ANTENNA_max);
        theta_max = generateRandomDouble(30, 80);
        parameters.push_back(radius);
        parameters.push_back(theta_max);
    }

    void generateCustomParameters(const std::vector<double>& _params) override {
        name = "spherical_cap";
        parameters.clear();
        radius = _params[0];
        theta_max = _params[1];
        parameters.push_back(radius);
        parameters.push_back(theta_max);
    }

    void generateBodyVerticies() override {
        vertices.clear();
        center = { 0, 0, 0 };
        double theta_max_rad = theta_max * std::numbers::pi / 180.0;
        int RINGS = std::max(4, static_cast<int>(radius * theta_max_rad / polygon_size));
        int SECTORS = std::max(4, static_cast<int>(2 * std::numbers::pi * radius * std::sin(theta_max_rad) / polygon_size));
        RINGS = std::min(RINGS, 100);
        SECTORS = std::min(SECTORS, 100);
        for (int i = 0; i <= RINGS; ++i) {
            double theta = i * theta_max_rad / RINGS;
            for (int j = 0; j <= SECTORS; ++j) {
                double phi = j * 2 * std::numbers::pi / SECTORS;
                vertices.push_back({
                    center.x + radius * std::sin(theta) * std::cos(phi),
                    center.y + radius * std::sin(theta) * std::sin(phi),
                    center.z + radius * std::cos(theta)
                    });
            }
        }
    }

    void generateComponentVerticies(const std::vector<Point>& _limits, int _place) override {
        vertices.clear();
        double halfLength = radius;
        switch (_place) {
        case 0: center = { _limits[0].x + halfLength + OFFSET, 0, 0 }; break;
        case 1: center = { 0, _limits[1].y - halfLength - OFFSET, 0 }; break;
        case 2: center = { _limits[2].x - halfLength - OFFSET, 0, 0 }; break;
        case 3: center = { 0, _limits[3].y + halfLength + OFFSET, 0 }; break;
        }
        double theta_max_rad = theta_max * std::numbers::pi / 180.0;
        int RINGS = std::max(4, static_cast<int>(radius * theta_max_rad / polygon_size));
        int SECTORS = std::max(4, static_cast<int>(2 * std::numbers::pi * radius * std::sin(theta_max_rad) / polygon_size));
        RINGS = std::min(RINGS, 100);
        SECTORS = std::min(SECTORS, 100);
        for (int i = 0; i <= RINGS; ++i) {
            double theta = i * theta_max_rad / RINGS;
            for (int j = 0; j <= SECTORS; ++j) {
                double phi = j * 2 * std::numbers::pi / SECTORS;
                vertices.push_back({
                    center.x + radius * std::sin(theta) * std::cos(phi),
                    center.y + radius * std::sin(theta) * std::sin(phi),
                    center.z + radius * std::cos(theta)
                    });
            }
        }
    }

    std::vector<Triangle> getTriangles() const override {
        std::vector<Triangle> triangles;
        double theta_max_rad = theta_max * std::numbers::pi / 180.0;
        int RINGS = std::max(4, static_cast<int>(radius * theta_max_rad / polygon_size));
        int SECTORS = std::max(4, static_cast<int>(2 * std::numbers::pi * radius * std::sin(theta_max_rad) / polygon_size));
        RINGS = std::min(RINGS, 100);
        SECTORS = std::min(SECTORS, 100);
        
        // Handle apex (i=0) specially - all vertices at ring 0 are the same point
        // Create triangles from apex to first ring
        for (int j = 0; j < SECTORS; ++j) {
            int apex = 0;  // All points in ring 0 are at the apex
            int idx3 = (SECTORS + 1) + j;
            int idx4 = idx3 + 1;
            triangles.push_back({ vertices[apex], vertices[idx3], vertices[idx4] });
            triangles.push_back({ vertices[apex], vertices[idx4], vertices[idx3] });
        }
        
        // Handle remaining rings normally
        for (int i = 1; i < RINGS; ++i) {
            for (int j = 0; j < SECTORS; ++j) {
                int idx1 = i * (SECTORS + 1) + j;
                int idx2 = idx1 + 1;
                int idx3 = (i + 1) * (SECTORS + 1) + j;
                int idx4 = idx3 + 1;
                triangles.push_back({ vertices[idx1], vertices[idx2], vertices[idx3] });
                triangles.push_back({ vertices[idx2], vertices[idx4], vertices[idx3] });
                triangles.push_back({ vertices[idx1], vertices[idx3], vertices[idx2] });
                triangles.push_back({ vertices[idx2], vertices[idx3], vertices[idx4] });
            }
        }
        return triangles;
    }

    std::vector<Point> computeTriangleNormals(const Triangle& triangle) const override {
        Point vec1 = { triangle.v2.x - triangle.v1.x, triangle.v2.y - triangle.v1.y, triangle.v2.z - triangle.v1.z };
        Point vec2 = { triangle.v3.x - triangle.v1.x, triangle.v3.y - triangle.v1.y, triangle.v3.z - triangle.v1.z };
        Point normal = {
            vec1.y * vec2.z - vec1.z * vec2.y,
            vec1.z * vec2.x - vec1.x * vec2.z,
            vec1.x * vec2.y - vec1.y * vec2.x
        };
        double length = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
        if (length > 0) {
            normal.x /= length;
            normal.y /= length;
            normal.z /= length;
        }
        return { normal, {-normal.x, -normal.y, -normal.z} };
    }

    void storeParameters(std::ofstream& out) const override {
        out << "\tRadius: " << radius << "\n";
        out << "\tTheta Max: " << theta_max << "\n";
    }
};

class Parallelepiped : public Shape {
    double width, height, depth;
public:
    void generateRandomParameters() override {
        name = "parallelepiped";
        parameters.clear();
        width = generateRandomDouble(BODY_min, BODY_max);
        height = generateRandomDouble(BODY_min, BODY_max);
        depth = generateRandomDouble(BODY_min, BODY_max);
        parameters.push_back(width);
        parameters.push_back(height);
        parameters.push_back(depth);
    }

    void generateRandomParametersSolar() override {
        name = "parallelepiped";
        parameters.clear();
        width = generateRandomDouble(SOLAR_min, SOLAR_max);
        height = generateRandomDouble(SOLAR_min, SOLAR_max);
        depth = generateRandomDouble(SOLAR_min, SOLAR_max);
        parameters.push_back(width);
        parameters.push_back(height);
        parameters.push_back(depth);
    }

    void generateCustomParameters(const std::vector<double>& _params) override {
        name = "parallelepiped";
        parameters.clear();
        width = _params[0];
        height = _params[1];
        depth = _params[2];
        parameters.push_back(width);
        parameters.push_back(height);
        parameters.push_back(depth);
    }

    const int getSolarClass() const override { return 0; }

    void generateBodyVerticies() override {
        vertices.clear();
        center = { 0, 0, 0 };
        double halfWidth = width / 2, halfHeight = height / 2, halfDepth = depth / 2;
        vertices = {
            {-halfWidth, halfHeight, halfDepth}, {halfWidth, halfHeight, halfDepth}, {halfWidth, -halfHeight, halfDepth}, {-halfWidth, -halfHeight, halfDepth},
            {-halfWidth, halfHeight, -halfDepth}, {halfWidth, halfHeight, -halfDepth}, {halfWidth, -halfHeight, -halfDepth}, {-halfWidth, -halfHeight, -halfDepth}
        };
    }

    void generateComponentVerticies(const std::vector<Point>& _limits, int _place) override {
        vertices.clear();
        double halfWidth = width / 2, halfHeight = height / 2, halfDepth = depth / 2;
        switch (_place) {
        case 0: center = { _limits[0].x + halfWidth + OFFSET, 0, 0 }; break;
        case 1: center = { 0, _limits[1].y - halfHeight - OFFSET, 0 }; break;
        case 2: center = { _limits[2].x - halfWidth - OFFSET, 0, 0 }; break;
        case 3: center = { 0, _limits[3].y + halfHeight + OFFSET, 0 }; break;
        }
        vertices = {
            {center.x - halfWidth, center.y + halfHeight, center.z + halfDepth},
            {center.x + halfWidth, center.y + halfHeight, center.z + halfDepth},
            {center.x + halfWidth, center.y - halfHeight, center.z + halfDepth},
            {center.x - halfWidth, center.y - halfHeight, center.z + halfDepth},
            {center.x - halfWidth, center.y + halfHeight, center.z - halfDepth},
            {center.x + halfWidth, center.y + halfHeight, center.z - halfDepth},
            {center.x + halfWidth, center.y - halfHeight, center.z - halfDepth},
            {center.x - halfWidth, center.y - halfHeight, center.z - halfDepth}
        };
    }

    std::vector<Triangle> getTriangles() const override {
        std::vector<Triangle> triangles;
        std::vector<std::vector<Point>> faces = {
            {vertices[0], vertices[3], vertices[2], vertices[1]}, // �������� ����� (Z = halfDepth)
            {vertices[4], vertices[5], vertices[6], vertices[7]}, // ������ ����� (Z = -halfDepth)
            {vertices[0], vertices[1], vertices[5], vertices[4]}, // ������� ����� (Y = halfHeight)
            {vertices[2], vertices[3], vertices[7], vertices[6]}, // ������ ����� (Y = -halfHeight)
            {vertices[1], vertices[5], vertices[6], vertices[2]}, // ������ ����� (X = halfWidth)
            {vertices[0], vertices[4], vertices[7], vertices[3]}  // ����� ����� (X = -halfWidth)
        };
        for (const auto& face : faces) {
            auto faceTriangles = triangulateFaceTriangles(face, polygon_size);
            triangles.insert(triangles.end(), faceTriangles.begin(), faceTriangles.end());
        }
        return triangles;
    }

    std::vector<Point> computeTriangleNormals(const Triangle& triangle) const override {
        Point vec1 = { triangle.v2.x - triangle.v1.x, triangle.v2.y - triangle.v1.y, triangle.v2.z - triangle.v1.z };
        Point vec2 = { triangle.v3.x - triangle.v1.x, triangle.v3.y - triangle.v1.y, triangle.v3.z - triangle.v1.z };
        Point normal = {
            vec1.y * vec2.z - vec1.z * vec2.y,
            vec1.z * vec2.x - vec1.x * vec2.z,
            vec1.x * vec2.y - vec1.y * vec2.x
        };
        double length = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
        if (length > 0) {
            normal.x /= length;
            normal.y /= length;
            normal.z /= length;
        }
        Point triangle_center = {
            (triangle.v1.x + triangle.v2.x + triangle.v3.x) / 3.0,
            (triangle.v1.y + triangle.v2.y + triangle.v3.y) / 3.0,
            (triangle.v1.z + triangle.v2.z + triangle.v3.z) / 3.0
        };
        Point vec_to_center = {
            triangle_center.x - center.x,
            triangle_center.y - center.y,
            triangle_center.z - center.z
        };
        double dot_product = normal.x * vec_to_center.x + normal.y * vec_to_center.y + normal.z * vec_to_center.z;
        if (dot_product < 0) {
            normal.x = -normal.x;
            normal.y = -normal.y;
            normal.z = -normal.z;
        }
        return { normal };
    }

    void storeParameters(std::ofstream& out) const override {
        out << "\tWidth: " << width << "\n";
        out << "\tHeight: " << height << "\n";
        out << "\tDepth: " << depth << "\n";
    }
};

class TrapezoidalPrism : public Shape {
    double top_width, bottom_width, height, depth;
private:
    // ��������������� ������� ��� �������������� ������ (�������� � ������)
    std::vector<Triangle> triangulateTrapezoidFace(const std::vector<Point>& face, double polygon_size, bool reverse) const {
        std::vector<Triangle> triangles;
        // ������� ������: {top_left, top_right, bottom_right, bottom_left}
        // ��� ��������: {0: (-1.5, 5, -3), 1: (1.5, 5, -3), 2: (3.5, -5, -3), 3: (-3.5, -5, -3)}
        // ��� ������: {4: (-1.5, 5, 3), 5: (1.5, 5, 3), 6: (3.5, -5, 3), 7: (-3.5, -5, 3)}
        double min_y = face[0].y, max_y = face[0].y;
        for (const auto& p : face) {
            min_y = std::min(min_y, p.y);
            max_y = std::max(max_y, p.y);
        }
        double dy = max_y - min_y;
        int ny = static_cast<int>(std::ceil(dy / polygon_size));
        int nx = static_cast<int>(std::ceil(std::max(std::abs(face[1].x - face[0].x), std::abs(face[2].x - face[3].x)) / polygon_size));
        if (nx < 1) nx = 1;
        if (ny < 1) ny = 1;

        // ������ �����
        std::vector<std::vector<Point>> grid(ny + 1, std::vector<Point>(nx + 1));
        for (int j = 0; j <= ny; ++j) {
            double t_y = static_cast<double>(j) / ny;
            double y = max_y - t_y * dy; // �� max_y (����) � min_y (���)
            // ������������ X: ����� ����� (top_left -> bottom_left), ������ ����� (top_right -> bottom_right)
            double x_left = face[0].x * (1 - t_y) + face[3].x * t_y;  // top_left (0) -> bottom_left (3)
            double x_right = face[1].x * (1 - t_y) + face[2].x * t_y; // top_right (1) -> bottom_right (2)
            for (int i = 0; i <= nx; ++i) {
                double t_x = static_cast<double>(i) / nx;
                double x = x_left + t_x * (x_right - x_left);
                grid[j][i] = { x, y, face[0].z }; // Z ����������
            }
        }

        // ��������� ������������
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                Point p1 = grid[j][i];
                Point p2 = grid[j][i + 1];
                Point p3 = grid[j + 1][i];
                Point p4 = grid[j + 1][i + 1];
                if (reverse) {
                    triangles.push_back({ p1, p3, p2 });
                    triangles.push_back({ p2, p3, p4 });
                }
                else {
                    triangles.push_back({ p1, p2, p3 });
                    triangles.push_back({ p2, p4, p3 });
                }
            }
        }
        return triangles;
    }

    // ��������������� ������� ��� ��������� ������ (������ � �����)
    std::vector<Triangle> triangulateSlantedFace(const std::vector<Point>& face, double polygon_size, bool reverse) const {
        std::vector<Triangle> triangles;
        int ny = static_cast<int>(std::ceil((face[0].y - face[2].y) / polygon_size)); // �� ������ (Y)
        int nz = static_cast<int>(std::ceil((face[2].z - face[1].z) / polygon_size)); // �� ������� (Z)
        if (ny < 1) ny = 1;
        if (nz < 1) nz = 1;

        std::vector<std::vector<Point>> grid(ny + 1, std::vector<Point>(nz + 1));
        for (int j = 0; j <= ny; ++j) {
            double t_y = static_cast<double>(j) / ny;
            double y = face[0].y * (1 - t_y) + face[2].y * t_y;
            for (int k = 0; k <= nz; ++k) {
                double t_z = static_cast<double>(k) / nz;
                double z = face[0].z * (1 - t_z) + face[2].z * t_z;
                double x0 = face[0].x * (1 - t_y) + face[2].x * t_y;
                double x1 = face[3].x * (1 - t_y) + face[1].x * t_y;
                double x = x0 * (1 - t_z) + x1 * t_z;
                grid[j][k] = { x, y, z };
            }
        }

        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                Point p1 = grid[j][k];
                Point p2 = grid[j][k + 1];
                Point p3 = grid[j + 1][k];
                Point p4 = grid[j + 1][k + 1];
                if (reverse) {
                    triangles.push_back({ p1, p3, p2 });
                    triangles.push_back({ p2, p3, p4 });
                }
                else {
                    triangles.push_back({ p1, p2, p3 });
                    triangles.push_back({ p2, p4, p3 });
                }
            }
        }
        return triangles;
    }

public:
    void generateRandomParameters() override {
        name = "trapezoidal_prism";
        parameters.clear();
        bottom_width = generateRandomDouble(BODY_min, BODY_max);
        top_width = generateRandomDouble(int(bottom_width / 2.), (int)bottom_width);
        height = generateRandomDouble(BODY_min, BODY_max);
        depth = generateRandomDouble(BODY_min, BODY_max);
        parameters.push_back(top_width);
        parameters.push_back(bottom_width);
        parameters.push_back(height);
        parameters.push_back(depth);
    }

    void generateCustomParameters(const std::vector<double>& _params) override {
        name = "trapezoidal_prism";
        parameters.clear();
        top_width = _params[0];
        bottom_width = _params[1];
        height = _params[2];
        depth = _params[3];
        parameters.push_back(top_width);
        parameters.push_back(bottom_width);
        parameters.push_back(height);
        parameters.push_back(depth);
    }

    void generateBodyVerticies() override {
        vertices.clear();
        center = { 0, 0, 0 };
        double halfHeight = height / 2, halfDepth = depth / 2;
        double halfTopWidth = top_width / 2, halfBottomWidth = bottom_width / 2;
        vertices = {
            {-halfTopWidth, halfHeight, -halfDepth},    // 0: ������� �������� �����
            {halfTopWidth, halfHeight, -halfDepth},     // 1: ������� �������� ������
            {halfBottomWidth, -halfHeight, -halfDepth}, // 2: ������ �������� ������
            {-halfBottomWidth, -halfHeight, -halfDepth},// 3: ������ �������� �����
            {-halfTopWidth, halfHeight, halfDepth},     // 4: ������� ������ �����
            {halfTopWidth, halfHeight, halfDepth},      // 5: ������� ������ ������
            {halfBottomWidth, -halfHeight, halfDepth},  // 6: ������ ������ ������
            {-halfBottomWidth, -halfHeight, halfDepth}  // 7: ������ ������ �����
        };
        
    }

    void generateComponentVerticies(const std::vector<Point>& _limits, int _place) override {
        vertices.clear();
        double halfHeight = height / 2, halfDepth = depth / 2;
        double halfTopWidth = top_width / 2, halfBottomWidth = bottom_width / 2;
        double halfLength = std::max(halfTopWidth, halfBottomWidth);
        switch (_place) {
        case 0: center = { _limits[0].x + halfLength + OFFSET, 0, 0 }; break;
        case 1: center = { 0, _limits[1].y - halfHeight - OFFSET, 0 }; break;
        case 2: center = { _limits[2].x - halfLength - OFFSET, 0, 0 }; break;
        case 3: center = { 0, _limits[3].y + halfHeight + OFFSET, 0 }; break;
        }
        vertices = {
            {center.x - halfTopWidth, center.y + halfHeight, center.z - halfDepth},
            {center.x + halfTopWidth, center.y + halfHeight, center.z - halfDepth},
            {center.x + halfBottomWidth, center.y - halfHeight, center.z - halfDepth},
            {center.x - halfBottomWidth, center.y - halfHeight, center.z - halfDepth},
            {center.x - halfTopWidth, center.y + halfHeight, center.z + halfDepth},
            {center.x + halfTopWidth, center.y + halfHeight, center.z + halfDepth},
            {center.x + halfBottomWidth, center.y - halfHeight, center.z + halfDepth},
            {center.x - halfBottomWidth, center.y - halfHeight, center.z + halfDepth}
        };
       
    }

    std::vector<Triangle> getTriangles() const override {
        std::vector<Triangle> triangles;
        std::vector<std::vector<Point>> faces = {
            {vertices[0], vertices[1], vertices[2], vertices[3]}, // �������� ����� (��������, Z = -halfDepth)
            {vertices[4], vertices[5], vertices[6], vertices[7]}, // ������ ����� (��������, Z = halfDepth)
            {vertices[0], vertices[4], vertices[5], vertices[1]}, // ������� ����� (�������������, Y = halfHeight)
            {vertices[3], vertices[2], vertices[6], vertices[7]}, // ������ ����� (�������������, Y = -halfHeight)
            {vertices[1], vertices[2], vertices[6], vertices[5]}, // ������ ����� (���������)
            {vertices[0], vertices[3], vertices[7], vertices[4]}  // ����� ����� (���������)
        };
        for (size_t i = 0; i < faces.size(); ++i) {
            bool isSlanted = (i == 4 || i == 5); // ������ � ����� �����
            bool isTrapezoid = (i == 0 || i == 1); // �������� � ������ �����
            std::vector<Triangle> faceTriangles;
            if (isTrapezoid) {
                faceTriangles = triangulateTrapezoidFace(faces[i], polygon_size, false);
            }
            else if (isSlanted) {
                faceTriangles = triangulateSlantedFace(faces[i], polygon_size, false);
            }
            else {
                faceTriangles = Shape::triangulateFaceTriangles(faces[i], polygon_size, false);
            }
            triangles.insert(triangles.end(), faceTriangles.begin(), faceTriangles.end());

        }
        return triangles;
    }

    std::vector<Triangle> triangulateFaceTriangles(const std::vector<Point>& face, double polygon_size, bool reverse = false) const override {
        // ��������� ��� �����
        bool isSlanted = (std::abs(face[0].x - face[2].x) > 1e-6 || std::abs(face[1].x - face[3].x) > 1e-6);
        bool isTrapezoid = (std::abs(face[0].z - face[1].z) < 1e-6 && std::abs(face[2].z - face[3].z) < 1e-6 &&
            std::abs(face[0].x - face[1].x) != std::abs(face[2].x - face[3].x));

        if (isTrapezoid) {
            return triangulateTrapezoidFace(face, polygon_size, reverse);
        }
        else if (isSlanted) {
            return triangulateSlantedFace(face, polygon_size, reverse);
        }
        else {
            return Shape::triangulateFaceTriangles(face, polygon_size, reverse);
        }
    }

    std::vector<Point> computeTriangleNormals(const Triangle& triangle) const override {
        Point vec1 = { triangle.v2.x - triangle.v1.x, triangle.v2.y - triangle.v1.y, triangle.v2.z - triangle.v1.z };
        Point vec2 = { triangle.v3.x - triangle.v1.x, triangle.v3.y - triangle.v1.y, triangle.v3.z - triangle.v1.z };
        Point normal = {
            vec1.y * vec2.z - vec1.z * vec2.y,
            vec1.z * vec2.x - vec1.x * vec2.z,
            vec1.x * vec2.y - vec1.y * vec2.x
        };
        double length = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
        if (length > 0) {
            normal.x /= length;
            normal.y /= length;
            normal.z /= length;
        }
        Point triangle_center = {
            (triangle.v1.x + triangle.v2.x + triangle.v3.x) / 3.0,
            (triangle.v1.y + triangle.v2.y + triangle.v3.y) / 3.0,
            (triangle.v1.z + triangle.v2.z + triangle.v3.z) / 3.0
        };
        Point vec_to_center = {
            triangle_center.x - center.x,
            triangle_center.y - center.y,
            triangle_center.z - center.z
        };
        double dot_product = normal.x * vec_to_center.x + normal.y * vec_to_center.y + normal.z * vec_to_center.z;
        if (dot_product < 0) {
            normal.x = -normal.x;
            normal.y = -normal.y;
            normal.z = -normal.z;
        }

        return { normal };
    }

    void storeParameters(std::ofstream& out) const override {
        out << "\tTop Width: " << top_width << "\n";
        out << "\tBottom Width: " << bottom_width << "\n";
        out << "\tHeight: " << height << "\n";
        out << "\tDepth: " << depth << "\n";
    }
};
