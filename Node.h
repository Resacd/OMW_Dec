#pragma once
#include <iostream>
#include "VelocitySet.h"
#include "include.h"
#include "Struct.h"
#include "Macro.h"


class Node
{



public:
	std::vector<double> m_distributions;
	std::vector<double> m_newDistributions;
	std::vector<double> m_eqdistributions;
	std::shared_ptr<Node> m_neighbours[_nLatNodes];
	int x_position;
	int y_position;
	int z_position;

	NodeType NodeType_;
	BoundaryType BoundaryType_;
        double density;
        std::vector<double>velocity;

	std::vector<double> velocity_gradient; 
	double wall_shear_stress;
	std::vector<double> accumulated_shear_stress;



	double accumulated_shear_stress_magnitude;
	double osi;
	double tawss;
        std::vector<double> normal;

	Node() = default;
	Node(int coord_x, int coord_y, int coord_z, NodeType nt);
	~Node();

	double Equilibrium(std::shared_ptr<VelocitySet> velSet, int dir);
	void Stream(int dir);
	void Initialize(std::shared_ptr<VelocitySet> velSet, int dir);
	virtual double Density();
	virtual std::vector<double> Velocity(std::shared_ptr<VelocitySet> velSet);
	virtual void MacroscopicValue(std::shared_ptr<VelocitySet> velSet);
	double Pressure(std::shared_ptr<VelocitySet> velSet);
	void BounceBack(int dir, int opp_dir);
	virtual void ApplyEqBc();
	virtual void ApplyBc(double &vel);


	struct Hash {
		std::size_t operator()(const std::shared_ptr<Node>& node) const {
			return std::hash<int>()(node->x_position) ^
				(std::hash<int>()(node->y_position) << 1) ^
				(std::hash<int>()(node->z_position) << 2);
		}
	};


	struct Equal {
		bool operator()(const std::shared_ptr<Node>& lhs, const std::shared_ptr<Node>& rhs) const {
			return lhs->x_position == rhs->x_position &&
				lhs->y_position == rhs->y_position &&
				lhs->z_position == rhs->z_position;
		}
	};
};



