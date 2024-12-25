#pragma once
#ifndef _READER_H_
#define _READER_H_
#include "Struct.h"
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <cmath>

class Reader {

public:

	Reader();
	~Reader();
	std::vector<Point> readOffFile(const std::string& filename);
	std::vector<Triangle> readSTLASCII(const std::string& filename);
	std::vector<Triangle> readSTLBinary(std::string filename);

};

#endif // !_READER_H_

