#pragma once

#include "VelocitySet.h"
#include "Subdomain.h"
#include "OutputFile.h"
#include "Boundary.h"
#include <vector>
#include "pulsatile.h"
#include <queue>



class Simulation
{

public:
	OutputFile file1;
	int myrank_, num_procs_;
	SubDomain SubDomain_;
	Pulsatile Pulsatile_;
	double pressure_;
	double m_relaxation;
        double m_omega;
	double m_dt;
        double total_runtime;
        double Cs;
	std::vector<double> wall_shear_stress;
	double InletVelocity_;
	Simulation(SubDomain& SubDomain, double tau, double pressure, double del_t);
	~Simulation();
	void InitializeShearStressAccumulators();
	void Collision();
	void Streaming();
	void ApplyEqBc();
	void ApplyZouHeBc();
	void ApplyOutletBc();
	void ApplyWallBc();
	void log_data(const std::string& prefix, int source, int destination, int iteration, const std::vector<double>& data);
	void ComputeVelocityGradientAtWall();
	void ComputeWallShearStress(int iteration);
	void ComputeOSIAndTAWSS(int iteration);
	double FindClosestLatticePoint(const Point& targetPoint, InterpolationAttribute attribute);
	void InterpolateAttributesToSTLVertices(std::vector<Triangle>& stlTriangles, InterpolationAttribute attribute);

	void writeDataToFile(const std::vector<double>& data, const std::string& baseFilename, int iteration);
	void FaceCommunicate(int iter);
	void EdgeCommunicate();
	double getInletVelocity(double x, double period);
	void Iterate();
	void SetRelaxation(double tau);
	void WriteOutput(int iter);

};


