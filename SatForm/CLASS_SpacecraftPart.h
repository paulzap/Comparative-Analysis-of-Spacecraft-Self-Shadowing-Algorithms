#pragma once
#include "CLASS_Shape.h"
#include <memory>
#include <vector>
#include <fstream>

class SpacecraftPart {
protected:
    std::unique_ptr<Shape> shape;

public:
    virtual void initialize(double _pol_size) = 0;
    virtual void generateVerticies(const std::vector<Point>& points = {}, int place = -1) = 0;
    virtual void generateRandomParameters() = 0;
    virtual void storeToFile(std::ofstream& out) const = 0;
    virtual const Shape& getShape() const { return *shape; }
    virtual ~SpacecraftPart() = default;

    std::vector<Triangle> getTriangles() const {
        return shape->getTriangles();
    }
};

class Body : public SpacecraftPart {
public:
    void initialize(double _pol_size) override {
        int randValue = rand() % 100;
        if (randValue < 25) shape = std::make_unique<Sphere>();
        else if (randValue < 50) shape = std::make_unique<Cube>();
        else if (randValue < 75) shape = std::make_unique<TrapezoidalPrism>();
        else shape = std::make_unique<Parallelepiped>();
        shape->setPolygonSize(_pol_size);
        generateRandomParameters();
    }

    void generateRandomParameters() override {
        shape->generateRandomParameters();
    }

    void generateVerticies(const std::vector<Point>& points = {}, int place = -1) override {
        shape->generateBodyVerticies();
    }

    void storeToFile(std::ofstream& out) const override {
        out << "Body Component:\n";
        if (dynamic_cast<const TrapezoidalPrism*>(shape.get())) out << "\tType: TrapezoidalPrism\n";
        else if (dynamic_cast<const Parallelepiped*>(shape.get())) out << "\tType: Parallelepiped\n";
        else if (dynamic_cast<const Sphere*>(shape.get())) out << "\tType: Sphere\n";
        else if (dynamic_cast<const Cube*>(shape.get())) out << "\tType: Cube\n";
        else out << "\tType: Unknown\n";
        out << "\tParameters:\n";
        shape->storeParameters(out);
    }
};

class SolarPanel : public SpacecraftPart {
public:
    void initialize(double _pol_size) override {
        shape = std::make_unique<Plane>();
        shape->setPolygonSize(_pol_size);
        generateRandomParameters();
    }

    void initialize(int panel_class, const std::vector<double>& _params) {
        shape = std::make_unique<Plane>();
        shape->generateCustomParameters(_params);
    }

    void generateRandomParameters() override {
        shape->generateRandomParametersSolar();
    }

    void generateVerticies(const std::vector<Point>& limits = {}, int place = -1) override {
        shape->generateComponentVerticies(limits, place);
    }

    void storeToFile(std::ofstream& out) const override {
        out << "Solar Panel Component:\n";
        out << "\tType: Plane\n";
        out << "\tParameters:\n";
        shape->storeParameters(out);
    }
};

class Antenna : public SpacecraftPart {
public:
    void initialize(double _pol_size) override {
        int antenna_class = rand() % 3;
        switch (antenna_class) {
        case 0: shape = std::make_unique<Sphere>(); break;
        case 1: shape = std::make_unique<CurvedSurface>(); break;
        case 2: shape = std::make_unique<SphericalCap>(); break;
        }
        shape->setPolygonSize(_pol_size);
        generateRandomParameters();
    }

    void generateRandomParameters() override {
        shape->generateRandomParametersAntenna();
    }

    void generateVerticies(const std::vector<Point>& limits = {}, int place = -1) override {
        shape->generateComponentVerticies(limits, place);
    }

    void storeToFile(std::ofstream& out) const override {
        out << "Antenna Component:\n";
        if (dynamic_cast<const Sphere*>(shape.get())) out << "\tType: Sphere\n";
        else if (dynamic_cast<const CurvedSurface*>(shape.get())) out << "\tType: CurvedSurface\n";
        else if (dynamic_cast<const SphericalCap*>(shape.get())) out << "\tType: SphericalCap\n";
        out << "\tParameters:\n";
        shape->storeParameters(out);
    }
};