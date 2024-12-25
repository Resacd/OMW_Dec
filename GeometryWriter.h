#pragma once
#ifndef _WRITER_H_
#define _WRITER_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <array>
#include "Struct.h"
#include <unordered_map>
#include <memory> 

class Writer {

public:

    Writer();
    ~Writer();
    void writePointoffToFile(const std::vector<Point>& points, const std::string& filename);
    void writePointsToFile(const std::vector<Point>& pointsBox, const std::string& filename);
    void writeBoxToFile(const std::vector<std::vector<std::vector<Point>>>& points, const std::string& filename);
    void writePointsToFile(const std::vector<std::vector<std::vector<Point>>>& pointsBox, const std::string& filename);
    void writeLabeledPointsToFile(std::vector< std::vector< std::vector<std::shared_ptr<Point>>>> &pointsBox, const std::string& filename);
};	

#endif 
