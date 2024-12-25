//#pragma once
#include "Utility.h"
#include <algorithm>
#include <iostream>
#include <fstream>


Utility::Utility() {


}

Utility::~Utility() {


}

bool Utility::isPointInBoundary(Point& targetPoint, const std::vector<Point>& boundaryPoints) {

    for (const Point& boundaryPoint : boundaryPoints) {

        const double epsilon = 1e-6; 

        if (std::abs(targetPoint.x - boundaryPoint.x) < epsilon &&
            std::abs(targetPoint.y - boundaryPoint.y) < epsilon &&
            std::abs(targetPoint.z - boundaryPoint.z) < epsilon) {
            return true;
        }
    }


    return false;
}


bool Utility::AreTrianglesEqual(const Triangle& a, const Triangle& b) {
    return (
        a.x1 == b.x1 && a.y1 == b.y1 && a.z1 == b.z1 &&
        a.x2 == b.x2 && a.y2 == b.y2 && a.z2 == b.z2 &&
        a.x3 == b.x3 && a.y3 == b.y3 && a.z3 == b.z3 &&
        a.nx == b.nx && a.ny == b.ny && a.nz == b.nz
        );
}

void Utility::RemoveDuplicateTriangles(std::vector<Triangle>& triangles) {
    triangles.erase(std::unique(triangles.begin(), triangles.end(), AreTrianglesEqual), triangles.end());
}


