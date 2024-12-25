#pragma once

#include "VelocitySet.h"
#include "Macro.h"
class LatticeModel : public VelocitySet
{
public:

	LatticeModel();
	~LatticeModel();
	std::vector<std::pair<std::vector<int>, std::vector<double>>> vectorOfPairs;
	int OppositeDirection(int dir) const override;
	std::vector<std::pair<std::vector<int>, std::vector<double>>> Pair0fNormalUnkDis();
};


