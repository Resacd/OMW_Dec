#pragma once
#include "Node.h"

class Halo : public Node
{
private:
	std::vector<double> m_inletVel;
	std::shared_ptr<VelocitySet> velSet;
public:
	static const double PI;
public:
	Halo(int coord_x, int coord_y, int coord_z);
	~Halo();
};