#pragma once
#ifndef _SUBDOMAIN_H_
#define _SUBDOMAIN_H_

#include <vector>
#include "Node.h"
#include "HaloLayer.h"
#include "mpi.h"
#include <omp.h>
#include "typedef.h"
#include "Struct.h"
#include "Boundary.h"
#include "Utility.h"
#include "ComplexGeometry.h"


class SubDomain {

private:
	int myrank_, num_procs_;

	
public:

	Writer Writer;
	Utility Utility_;
	Geometry Geometry_;
	int start_x;
	int	start_y;
	int	start_z;
	int	end_x;
	int	end_y;
	int	end_z;
	int my_Nx;
	int my_Ny;
	int my_Nz;
	int halo_size;
	std::vector<int> gridSize_;
	std::vector<int> NewgridSize_;
	std::vector<int> m_subdomainSize;
	int SubdomainSize_;
	std::shared_ptr<VelocitySet> m_velSet;
	std::vector< std::vector< std::vector<std::shared_ptr<Node>>>> m_lattice;
	std::vector<std::shared_ptr<Node>> HaloLattice_;
	std::vector<std::shared_ptr<Node>> Involved_Lat_;
	std::vector<std::shared_ptr<Node>> BoundaryNode_;
	std::vector<std::shared_ptr<Node>> InletBoundaryNode_;
	std::vector<std::shared_ptr<Node>> OutletBoundaryNode1_;
	std::vector<std::shared_ptr<Node>> OutletBoundaryNode2_;
	std::vector<std::shared_ptr<Node>> WallBoundaryNode_;
	std::vector<std::shared_ptr<Node>> nextToTheWall;

	int my_up, my_down, my_left, my_right, my_front, my_back;
	int right_up_, right_down_, left_up_, left_down_, up_front_, up_back_;
	int down_front_, down_back_, right_front_, right_back_, left_front_, left_back_;
	std::vector<Neighbors> neighbors;
	MPI_Comm cart_comm;
	std::vector<int> dims;
	std::vector<int> coords;
	std::vector<int>need_communication;

	std::vector<std::pair<std::tuple<int, int, int>, std::tuple<int, int, int>>> local2global_;
	SubDomain() = default;
	SubDomain(Geometry& Geometry);
	~SubDomain();
	void InitSubDomain();
	void CreateNodes();
	bool IsInSubDomain(int coord_x, int coord_y, int coord_z);
	void Connectivity();
	void IdentifyNextToWallLattices();
        double CalculateAngle(const Triangle& triangle, double x, double y, double z);

        double CalculateVertexAngle(double vx, double vy, double vz, double v1x, double v1y, double v1z,
                                double v2x, double v2y, double v2z);


	void calculateBoundaryNormals();
	void PrintNonNullNeighbors();
	void Initialize();
	void Need_Comm();
	void GetFluidHaloLat();
	const std::shared_ptr<VelocitySet> GetVelSet();
	void WriteHaloLatticeToVTK();
	void WriteFluidLatticeToVTK();
	void WriteLatticeToVTK(const std::string& filename, const std::vector<std::shared_ptr<Node>>& lattice);

	friend class Simulation;
	friend class Node;


struct NodeGlobalHasher {
    Geometry& geometry_;

    NodeGlobalHasher(Geometry& geometry) : geometry_(geometry) {}

    size_t operator()(const std::shared_ptr<Node>& node) const {
        // Use the global index for hashing
        std::vector<int> globalIndex = geometry_.Local2Global(node->x_position, node->y_position, node->z_position);
        return std::hash<int>()(globalIndex[0]) ^ std::hash<int>()(globalIndex[1]) ^ std::hash<int>()(globalIndex[2]);
    }
};

struct NodeGlobalEqual {
    Geometry& geometry_;

    NodeGlobalEqual(Geometry& geometry) : geometry_(geometry) {}

    bool operator()(const std::shared_ptr<Node>& lhs, const std::shared_ptr<Node>& rhs) const {
        // Compare the global indices for equality
        std::vector<int> globalIndexLhs = geometry_.Local2Global(lhs->x_position, lhs->y_position, lhs->z_position);
        std::vector<int> globalIndexRhs = geometry_.Local2Global(rhs->x_position, rhs->y_position, rhs->z_position);
        return globalIndexLhs == globalIndexRhs;
    }
};




};
#endif


