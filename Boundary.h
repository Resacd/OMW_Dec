#pragma once
#include "Node.h"
#include "Macro.h"
#include "Utility.h"

class Boundary : public Node
{

private:
	std::vector<double> m_inletVel;
	std::shared_ptr<VelocitySet> velSet;


public:
	Utility Utility;
	Boundary(int coord_x, int coord_y, int coord_z, std::shared_ptr <VelocitySet> velset,NodeType nt);
	~Boundary();
	void SetInletVel(double ux, double uy,double uz);
	double Density() override;
	virtual std::vector<double> Velocity(std::shared_ptr<VelocitySet> velSet) override;
	void ApplyEqBc() override;
	void ApplyBc(double& vel) override;
	void MacroscopicValue(std::shared_ptr<VelocitySet> velSet) override;

	

};

