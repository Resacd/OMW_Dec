#pragma once
#ifndef _UTILITY_H_
#define _UTILITY_H_

#include <iostream>
#include <vector>
#include <unordered_set>
#include "Struct.h"
#include <tuple>
#include <map>
#include <set>
#include <fstream>

class Utility {

public:
    Utility();
    ~Utility();
    bool isPointInBoundary(Point& targetPoint, const std::vector<Point>& boundaryPoints);
   static bool AreTrianglesEqual(const Triangle& a, const Triangle& b);
    void RemoveDuplicateTriangles(std::vector<Triangle>& triangles);
    void onMouse(int event, int x, int y, int flags, void* userdata);
    void processImage(const std::string& imageFile);


};



#endif // !_UTILITY_
