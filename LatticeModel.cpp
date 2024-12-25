#ifndef LATTICE_MODEL_H
#define LATTICE_MODEL_H

#include "LatticeModel.h"
#include <iostream>



#ifdef D2Q9

LatticeModel::LatticeModel()
	: VelocitySet(
		{ 4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0 , 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0 },
		{ {0, 0}, {1, 0}, {0, 1}, {-1, 0}, {0, -1}, {1, 1}, {-1, 1}, {-1, -1}, {1, -1} },
		1.0 / 3.0
	)
{
}

int LatticeModel::OppositeDirection(int dir) const
{
	switch (dir)
	{
	case 0:
		return 0;
	case 1:
		return 3;
	case 2:
		return 4;
	case 3:
		return 1;
	case 4:
		return 2;
	case 5:
		return 7;
	case 6:
		return 8;
	case 7:
		return 5;
	case 8:
		return 6;
	}
}

#endif



#ifdef D3Q19

LatticeModel::LatticeModel()
	: VelocitySet(
		{ 2.0 / 36.0,      //0
		  2.0 / 36.0,      //1 
		  2.0 / 36.0,      //2
		  2.0 / 36.0,      //3
		  2.0 / 36.0,      //4
		  2.0 / 36.0,      //5
		  1.0 / 36.0,      //6
		  1.0 / 36.0,      //7
		  1.0 / 36.0,      //8
		  1.0 / 36.0,      //9
		  1.0 / 36.0,      //10
		  1.0 / 36.0,      //11
		  1.0 / 36.0,      //12
		  1.0 / 36.0,      //13
		  1.0 / 36.0,      //14
		  1.0 / 36.0,      //15
		  1.0 / 36.0,      //16
		  1.0 / 36.0,      //17
		  12.0 / 36.0,     //18

		},
		{ { 1,  0, 0},     //0
		  {-1,  0, 0},     //1
		  { 0,  1, 0},     //2
		  { 0, -1, 0},     //3
		  { 0,  0, 1},     //4 
		  { 0,  0,-1},     //5 
		  { 1,  1, 0},     //6
		  { 1, -1, 0},     //7   
		  { 1,  0, 1},     //8
		  { 1,  0,-1},     //9
		  {-1,  1, 0},     //10
		  {-1, -1, 0},     //11
		  {-1,  0, 1},     //12
		  {-1,  0,-1},     //13
		  { 0,  1, 1},     //14
		  { 0,  1,-1},     //15
		  { 0, -1, 1},     //16
		  { 0, -1,-1},     //17
		  { 0,  0, 0},     //18
		},


		{ { {0, 6, 7, 8, 9}, {1, 0, 0} },
		{ {1, 10, 11, 12, 13}, {-1, 0, 0} },
		{ {3, 7, 11, 16, 17}, {0, -1, 0} },
		{ {2, 6, 10, 14, 15}, {0, 1, 0} },
		{ {4, 8, 12, 14, 16}, {0, 0, 1} },
		{ {5, 9, 13, 15, 17}, {0, 0, -1} } },

		1.0 / 3.0
	)
{
}

int LatticeModel::OppositeDirection(int dir) const
{
	switch (dir)
	{
	case 0:
		return 1;
	case 1:
		return 0;
	case 2:
		return 3;
	case 3:
		return 2;
	case 4:
		return 5;
	case 5:
		return 4;
	case 6:
		return 11;
	case 7:
		return 10;
	case 8:
		return 13;
	case 9:
		return 12;
	case 10:
		return 7;
	case 11:
		return 6;
	case 12:
		return 9;
	case 13:
		return 8;
	case 14:
		return 17;
	case 15:
		return 16;
	case 16:
		return 15;
	case 17:
		return 14;
	case 18:
		return 18;
	}
}


std::vector<std::pair<std::vector<int>, std::vector<double>>> LatticeModel::Pair0fNormalUnkDis() {


	std::vector<int> EastSide{ 0,6,7,8,9 };
	std::vector<double> EastNormal{ 1,0,0 };
	std::pair<std::vector<int>, std::vector<double>> EastPair = std::make_pair(EastSide, EastNormal);

	std::vector<int> WestSide{ 1,10,11,12,13 };
	std::vector<double> WestNormal{ -1,0,0 };
	std::pair<std::vector<int>, std::vector<double>> WestPair = std::make_pair(WestSide, WestNormal);

	std::vector<int> FrontSide{ 3,7,11,16,17 };
	std::vector<double> FrontNormal{ 0,-1,0 };
	std::pair<std::vector<int>, std::vector<double>> FrontPair = std::make_pair(FrontSide, FrontNormal);

	std::vector<int> BackSide{ 2,6,10,14,15 };
	std::vector<double> BackNormal{ 0, 1,0 };
	std::pair<std::vector<int>, std::vector<double>> BackPair = std::make_pair(BackSide, BackNormal);

	std::vector<int> UpSide{ 4,8,12,14,16 };
	std::vector<double> UpNormal{ 0, 0,1 };
	std::pair<std::vector<int>, std::vector<double>> UpPair = std::make_pair(UpSide, UpNormal);

	std::vector<int> DownSide{ 5,9,13,15,17 };
	std::vector<double> DownNormal{ 0, 0,-1 };
	std::pair<std::vector<int>, std::vector<double>> DownPair = std::make_pair(DownSide, DownNormal);

	vectorOfPairs.push_back(EastPair);
	vectorOfPairs.push_back(WestPair);
	vectorOfPairs.push_back(FrontPair);
	vectorOfPairs.push_back(BackPair);
	vectorOfPairs.push_back(UpPair);
	vectorOfPairs.push_back(DownPair);

	for (const auto& pair : vectorOfPairs) {
		std::cout << "Pair:\n";
		std::cout << "First vector: ";
		for (const auto& num : pair.first) {
			std::cout << num << " ";
		}
		std::cout << "\nSecond vector: ";
		for (const auto& num : pair.second) {
			std::cout << num << " ";
		}
		std::cout << "\n\n";
	}

	return vectorOfPairs;
}

#endif

#ifdef D3Q27

LatticeModel::LatticeModel()
	: VelocitySet(
		{ 8.0 / 27.0,      //0
		  2.0 / 27.0,      //1 
		  2.0 / 27.0,      //2
		  2.0 / 27.0,      //3
		  2.0 / 27.0,      //4
		  2.0 / 27.0,      //5
		  2.0 / 27.0,      //6
		  1.0 / 54.0,      //7
		  1.0 / 54.0,      //8
		  1.0 / 54.0,      //9
		  1.0 / 54.0,      //10
		  1.0 / 54.0,      //11
		  1.0 / 54.0,      //12
		  1.0 / 54.0,      //13
		  1.0 / 54.0,      //14
		  1.0 / 54.0,      //15
		  1.0 / 54.0,      //16
		  1.0 / 54.0,      //17
		  1.0 / 54.0,      //18
		  1.0 / 216.0,     //19
		  1.0 / 216.0,     //20
		  1.0 / 216.0,     //21
		  1.0 / 216.0,     //22
		  1.0 / 216.0,     //23
		  1.0 / 216.0,     //24
		  1.0 / 216.0,     //25
		  1.0 / 216.0,     //26

		},
		{ { 0,  0, 0},     //0
		  { 1,  0, 0},     //1
		  {-1,  0, 0},     //2
		  { 0,  1, 0},     //3
		  { 0, -1, 0},     //4 
		  { 0,  0, 1},     //5 
		  { 0,  0,-1},     //6
		  { 1,  1, 0},     //7   
		  { 1, -1, 0},     //8
		  {-1,  1, 0},     //9
		  {-1, -1, 0},     //10
		  { 1,  0, 1},     //11
		  {-1,  0, 1},     //12
		  { 1,  0,-1},     //13
		  {-1,  0,-1},     //14
		  { 0,  1, 1},     //15
		  { 0,  1,-1},     //16
		  { 0, -1, 1},     //17
		  { 0, -1,-1},     //18
		  { 1,  1, 1},     //19
		  {-1, -1,-1},     //20
		  { 1,  1,-1},     //21
		  {-1, -1, 1},     //22
		  { 1, -1, 1},     //23
		  {-1,  1,-1},     //24
		  { 1, -1,-1},     //25
		  {-1,  1, 1},     //26


		},
		1.0 / 3.0
	)
{
}

int LatticeModel::OppositeDirection(int dir) const
{
	switch (dir)
	{
	case 0:
		return 0;
	case 1:
		return 2;
	case 2:
		return 1;
	case 3:
		return 4;
	case 4:
		return 3;
	case 5:
		return 6;
	case 6:
		return 5;
	case 7:
		return 10;
	case 8:
		return 9;
	case 9:
		return 8;
	case 10:
		return 7;
	case 11:
		return 14;
	case 12:
		return 13;
	case 13:
		return 12;
	case 14:
		return 11;
	case 15:
		return 18;
	case 16:
		return 17;
	case 17:
		return 16;
	case 18:
		return 15;
	case 19:
		return 20;
	case 20:
		return 19;
	case 21:
		return 22;
	case 22:
		return 21;
	case 23:
		return 24;
	case 24:
		return 23;
	case 25:
		return 26;
	case 26:
		return 25;

	}
}


#endif

	LatticeModel::~LatticeModel()
	{
	}

#endif

