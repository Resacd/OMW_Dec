#pragma once
#ifndef _MESH_GENERATOR_H_
#define _MESH_GENERATOR_H_

#include "Struct.h"
#include <tuple>
#include <vector>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <fstream>
#include <mpi.h>
#include <omp.h>
#include <map>
#include "GeometryWriter.h"
#include "GeometryReader.h"
#include "mpi_call.h"
#include "Utility.h"
#include "VelocitySet.h"
#include "Macro.h"
#include <filesystem>
# include <memory>
#include <sys/stat.h> 



class Geometry {

public:
    Writer Writer;
    Reader Reader;
    Utility Utility;
    MPI_CALL<MPI_Datatype> mpi_caller;
    std::vector<std::vector<std::vector<std::shared_ptr<Point>>>> pointsBox;
    std::vector<Triangle> globalTriangles;
    std::vector<Triangle> localTriangles;
    std::string stlFileName;
    BoundingBox boundingbox;
    BoundingBox subBox;
    Octree octree;
    int maxDepth;
    std::vector<Point> boundaryPoints;
    std::vector<Point> InletboundaryPoints;
    std::vector<Point> OutletboundaryPoints1;
    std::vector<Point> OutletboundaryPoints2;
    std::vector<Point> WallboundaryPoints;
    std::vector<Point> NextToTheWall;
    int globalNumtris;
    std::vector<int> GeomSize;
    double spacing;
    double rayLength;
    int myrank, num_procs;
    MPI_Comm cart_comm;
    int halo_size;

    int nx;
    int ny;
    int nz;


    int my_Nx;
    int my_Ny;
    int my_Nz;

    int local_end_x;
    int local_start_x;

    int local_end_y;
    int local_start_y;

    int local_end_z;
    int local_start_z;

    int dims[3], coords[3];
    int my_left, my_right, my_front, my_back, my_up, my_down;
    std::vector<Neighbors> neighbors;

    int right_up_, right_down_, left_up_, left_down_, up_front_, up_back_;
    int down_front_, down_back_, right_front_, right_back_, left_front_, left_back_;

    std::shared_ptr<VelocitySet> m_velSet;
    char raydirection;

    Geometry() = default;
    Geometry(const std::string& FileName, std::shared_ptr<VelocitySet> velSet, double spacing, char raydirection);
    ~Geometry();
    void GetSubBoundingBox();
    void Create_cart_comm();
    void getNeighborInfo(int neighbor_rank, Neighbors& neighbor);
    bool is_valid_coordinate(std::vector<int>& coords);
    int get_rank(MPI_Comm comm, const std::vector<int>& coords);
    void Add_Ghost();
    bool DoesTriangleIntersectBox(const Triangle& tri, const BoundingBox& box);
    void DistributeTriangles();
    bool IsGhost(int i, int j, int k);
    void generatePointsInSubBox();
    bool IsInBox(int coord_x, int coord_y, int coord_z);
    std::vector<int>RetrieveLocalInd(int i, int j, int k);
    std::vector<int>Local2Global(int i, int j, int k);
    bool doLinesIntersect(const Line& line1, const Line& line2);
    bool isPointInsideTriangle(double x, double y, const Triangle& triangle);
    int countIntersectinsTriangle(Line& line, std::vector<Triangle>& posintersectedtri);
    void adjustLineCoordinates(Line& line, BoundingBox& boundingbox, double rayOffset, int directionSwitch);
    void labelPoints();
    void modifier();
    void Comm_Ghost();
    void OutlierHandler();
    std::vector<Point> findBoundaryPoints();
    void classifyBoundaryPoints();
    void writePointsBox(const std::string& filename);
    void writeLabeledPointsToFile(const std::string& filename);
    void writePointsToFile(const std::vector<Point>& points, const std::string& foldername, const std::string& filename);
};

#endif


