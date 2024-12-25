#pragma once

#ifndef _USER_DEFINED_STRUCTS_H_
#define _USER_DEFINED_STRUCTS_H_
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <map>
#include <string>

enum class PointBoxType {
    real,
    ghost

};



enum class BoundaryType {
    Inlet,
    Outlet,
    Wall,
    None
};

enum class NodeType {

    Fluid,
    Boundary,
    NotFluid,
    Halo
};

enum class InterpolationAttribute {
    WSS,
    OSI,
    TAWSS
};


inline std::string nodeTypeToString(NodeType type) {
    switch (type) {
    case NodeType::Fluid: return "Fluid";
    case NodeType::Boundary: return "InBoundary";
    case NodeType::NotFluid: return "NotFluid";
    case NodeType::Halo: return "Halo";
    default: return "Unknown";
    }
}



struct Point {

    double x;
    double y;
    double z;

    PointBoxType PointBoxType;
    NodeType NodeType;


    double shearStress;

    Point operator-(const Point& other) const {
        return { x - other.x, y - other.y, z - other.z };
    }

    Point operator+(const Point& other) const {
        return { x + other.x, y + other.y, z + other.z };
    }

    Point operator*(const double scalar) const {
        return { x * scalar, y * scalar, z * scalar };
    }

    bool operator==(const Point& other) const {
        return x == other.x && y == other.y && z == other.z;
    }

    bool operator<(const Point& other) const {
        if (x != other.x) return x < other.x;
        if (y != other.y) return y < other.y;
        return z < other.z;
    }
};



struct PointComparator {
    bool operator()(const Point& lhs, const Point& rhs) const {
        if (lhs.x != rhs.x) return lhs.x < rhs.x;
        if (lhs.y != rhs.y) return lhs.y < rhs.y;
        return lhs.z < rhs.z;
    }
};

struct Line {
    Point start;
    Point end;

    Point direction() const {
        return end - start;
    }
};


struct Triangle {
    double x1, y1, z1;
    double x2, y2, z2;
    double x3, y3, z3;
    double nx, ny, nz;
    double wss1, wss2, wss3;
    double tawss1, tawss2, tawss3;
    double osi1, osi2, osi3;
};

struct BoundingBox {
    double min_x;
    double max_x;
    double min_y;
    double max_y;
    double min_z;
    double max_z;
};

struct Edge {
    double x1, y1, z1;
    double x2, y2, z2;

    Edge(double x1, double y1, double z1, double x2, double y2, double z2)
        : x1(x1), y1(y1), z1(z1), x2(x2), y2(y2), z2(z2) {}

    bool operator==(const Edge& other) const {
        return (x1 == other.x1 && y1 == other.y1 && z1 == other.z1 && x2 == other.x2 && y2 == other.y2 && z2 == other.z2) ||
            (x1 == other.x2 && y1 == other.y2 && z1 == other.z2 && x2 == other.x1 && y2 == other.y1 && z2 == other.z1);
    }
};


namespace std {
    template <>
    struct hash<Edge> {
        std::size_t operator()(const Edge& e) const {
            std::size_t h1 = std::hash<double>()(e.x1);
            std::size_t h2 = std::hash<double>()(e.y1);
            std::size_t h3 = std::hash<double>()(e.z1);
            std::size_t h4 = std::hash<double>()(e.x2);
            std::size_t h5 = std::hash<double>()(e.y2);
            std::size_t h6 = std::hash<double>()(e.z2);
            return h1 ^ h2 ^ h3 ^ h4 ^ h5 ^ h6;
        }
    };
}

struct Octree {
    BoundingBox bBox;
    std::vector<Triangle> tri; 
    std::vector<Octree> children;    
};


struct NodeData {
    int global_x;
    int global_y;
    int global_z;
    double velX;
    double velY;
    double velZ;
    double velMag;
};

struct PairHash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ h2;
    }
};


struct Neighbors {
    int rank;
    int coords[3];
};

#endif // !_USER_DEFINED_STRUCTS_H_


