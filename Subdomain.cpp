#include "Subdomain.h"

SubDomain::SubDomain(Geometry& Geometry) {

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs_);
    Geometry_ = Geometry;
    m_velSet =Geometry_.m_velSet;
    InitSubDomain();
    CreateNodes();
    GetFluidHaloLat();
    //WriteFluidLatticeToVTK();
    //WriteHaloLatticeToVTK();
    Connectivity();
    IdentifyNextToWallLattices();
    calculateBoundaryNormals();
    Initialize();
    Need_Comm();
    MPI_Comm_free(&cart_comm);  
    //PrintNonNullNeighbors();

   // std::string wallFileName = "wallboundaryLattice_" + std::to_string(myrank_) + ".vtk";
   // std::string inletFileName = "inletboundaryLattice_" + std::to_string(myrank_) + ".vtk";
    //std::string outletFileName1 = "outletboundaryLattice1_" + std::to_string(myrank_) + ".vtk";
    //std::string outletFileName2 = "outletboundaryLattice2_" + std::to_string(myrank_) + ".vtk";

   // WriteLatticeToVTK(inletFileName.c_str(), InletBoundaryNode_);
    //WriteLatticeToVTK(outletFileName1.c_str(), OutletBoundaryNode_);
   // WriteLatticeToVTK(wallFileName.c_str(), WallBoundaryNode_);


}


SubDomain::~SubDomain() {

}


void SubDomain::InitSubDomain() {

    cart_comm = Geometry_.cart_comm;
    neighbors = Geometry_.neighbors;
    dims = { Geometry_.dims[0], Geometry_.dims[1], Geometry_.dims[2] };
    coords = { Geometry_.coords[0], Geometry_.coords[1], Geometry_.coords[2] };

    my_left = Geometry_.my_left;
    my_right = Geometry_.my_right;
    my_front = Geometry_.my_front;
    my_back = Geometry_.my_back;
    my_up = Geometry_.my_up;
    my_down = Geometry_.my_down;

    my_Nx = Geometry_.my_Nx;
    my_Ny = Geometry_.my_Ny;
    my_Nz = Geometry_.my_Nz;

    start_x = Geometry_.local_start_x;
    start_y = Geometry_.local_start_y;
    start_z = Geometry_.local_start_z;

    end_x = Geometry_.local_end_x;
    end_y = Geometry_.local_end_y;
    end_z = Geometry_.local_end_z;


    //for (int p = 0; p < num_procs_; p++) {
    //    if (myrank_ == p) {
    //        printf("Processor %d:\n", myrank_);
    //        printf("my_Nx = %d, my_Ny = %d, my_Nz = %d\n", my_Nx, my_Ny, my_Nz);
    //        printf("start_x = %d, end_x = %d\n", start_x, end_x);
    //        printf("start_y = %d, end_y = %d\n", start_y, end_y);
    //        printf("start_z = %d, end_z = %d\n", start_z, end_z);
    //        fflush(stdout);
    //    }
    //    MPI_Barrier(cart_comm_);
    //}
}

bool SubDomain::IsInSubDomain(int coord_x, int coord_y, int coord_z) {

    if (coord_x >= 0 && coord_x < my_Nx
        && coord_y >= 0 && coord_y < my_Ny
        && coord_z >= 0 && coord_z < my_Nz)
        return true;
    else
        return false;

}

void SubDomain::CreateNodes() {

    std::vector<Point> InletBoundary = Geometry_.InletboundaryPoints;
    std::vector<Point> OutletBoundary1 = Geometry_.OutletboundaryPoints1;
    std::vector<Point> OutletBoundary2 = Geometry_.OutletboundaryPoints2;
    std::vector<Point> WallBoundary = Geometry_.WallboundaryPoints;


    m_lattice.resize(my_Nx);

    for (int i = 0; i < my_Nx; i++) {
        m_lattice[i].resize(my_Ny);
        for (int j = 0; j < my_Ny; j++) {
            m_lattice[i][j].resize(my_Nz);
            for (int k = 0; k < my_Nz; k++) {
                m_lattice[i][j][k] = nullptr;
            }
        }
    }

    int local_wallboundary_count = 0;

    for (int i = 0; i < my_Nx; i++) {
        for (int j = 0; j < my_Ny; j++) {
            for (int k = 0; k < my_Nz; k++) {
                if (Geometry_.pointsBox[i][j][k]->PointBoxType == PointBoxType::ghost) {
                    m_lattice[i][j][k] = std::make_shared<Node>(i, j, k, NodeType::Halo);
                }
                else if (Geometry_.pointsBox[i][j][k]->PointBoxType == PointBoxType::real && Geometry_.pointsBox[i][j][k]->NodeType == NodeType::Fluid) {

                    std::shared_ptr<Point> GlobalcurrentPointPtr = Geometry_.pointsBox[i][j][k];
                    Point& GlobalcurrentPoint = *GlobalcurrentPointPtr;


                    if (Utility_.isPointInBoundary(GlobalcurrentPoint, InletBoundary)) {

                        m_lattice[i][j][k] = std::make_shared<Boundary>(i, j, k, m_velSet, NodeType::Boundary);
                        m_lattice[i][j][k]->BoundaryType_ = BoundaryType::Inlet;
                        InletBoundaryNode_.push_back(m_lattice[i][j][k]);
                    }

                    else if (Utility_.isPointInBoundary(GlobalcurrentPoint, OutletBoundary1)) {

                        m_lattice[i][j][k] = std::make_shared<Node>(i, j, k, NodeType::Boundary);
                        m_lattice[i][j][k]->BoundaryType_ = BoundaryType::Outlet;
                        OutletBoundaryNode1_.push_back(m_lattice[i][j][k]);
                    }

                    else if (Utility_.isPointInBoundary(GlobalcurrentPoint, OutletBoundary2)) {

                        m_lattice[i][j][k] = std::make_shared<Node>(i, j, k, NodeType::Boundary);
                        m_lattice[i][j][k]->BoundaryType_ = BoundaryType::Outlet;
                        OutletBoundaryNode2_.push_back(m_lattice[i][j][k]);
                    }


                    else if (Utility_.isPointInBoundary(GlobalcurrentPoint, WallBoundary)) {

                        m_lattice[i][j][k] = std::make_shared<Node>(i, j, k, NodeType::Boundary);
                        m_lattice[i][j][k]->BoundaryType_ = BoundaryType::Wall;
                        WallBoundaryNode_.push_back(m_lattice[i][j][k]);
                        local_wallboundary_count++;
                    }

                    else {
                        m_lattice[i][j][k] = std::make_shared<Node>(i, j, k, NodeType::Fluid);
                    }
                }
            }
        }
    }

    int total_wallboundary_count = 0;
    MPI_Reduce(&local_wallboundary_count, &total_wallboundary_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        std::cout << "Total number of boundary nodes: " << total_wallboundary_count << std::endl;
    }

}


void SubDomain::Connectivity() {

    for (int i = 0; i < my_Nx; i++) {
        for (int j = 0; j < my_Ny; j++) {
            for (int k = 0; k < my_Nz; k++) {
                if (m_lattice[i][j][k] != nullptr) {

                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        int x_neighbor = i + m_velSet->GetDirection(dir)[0];
                        int y_neighbor = j + m_velSet->GetDirection(dir)[1];
                        int z_neighbor = k + m_velSet->GetDirection(dir)[2];

                        if (IsInSubDomain(x_neighbor, y_neighbor, z_neighbor) &&
                            m_lattice[x_neighbor][y_neighbor][z_neighbor] != nullptr) {
                            m_lattice[i][j][k]->m_neighbours[dir] = m_lattice[x_neighbor][y_neighbor][z_neighbor];
                        }
                        else {
                            m_lattice[i][j][k]->m_neighbours[dir] = nullptr;
                        }
                    }
                }
            }
        }
    }
}




void SubDomain::IdentifyNextToWallLattices() {
    std::unordered_set<std::shared_ptr<Node>, NodeGlobalHasher, NodeGlobalEqual> uniqueNextToWall(0, NodeGlobalHasher(Geometry_), NodeGlobalEqual(Geometry_));


    for (auto& wallNode : WallBoundaryNode_) {
        for (int dir = 0; dir < _nLatNodes; dir++) {
            auto neighbor = wallNode->m_neighbours[dir];

            if (neighbor) {
                if (neighbor->NodeType_ == NodeType::Fluid && neighbor->BoundaryType_ != BoundaryType::Wall) {
                    uniqueNextToWall.insert(neighbor);
                }
            }
        }
    }

    nextToTheWall.assign(uniqueNextToWall.begin(), uniqueNextToWall.end());
}






/*
void SubDomain::IdentifyNextToWallLattices() {

    std::unordered_set<std::shared_ptr<Node>, Node::Hash, Node::Equal> uniqueNextToWall;


    for (auto& wallNode : WallBoundaryNode_) {
        for (int dir = 0; dir < _nLatNodes; dir++) {
            auto neighbor = wallNode->m_neighbours[dir];

            if (neighbor) {

                if (neighbor->NodeType_ == NodeType::Fluid && neighbor->BoundaryType_ != BoundaryType::Wall) {

                    uniqueNextToWall.insert(neighbor);
                }
            }
        }
    }

    nextToTheWall.assign(uniqueNextToWall.begin(), uniqueNextToWall.end());
}

*/

//Voronoi area weighting

void SubDomain::calculateBoundaryNormals() {
    for (auto& wall_lat : nextToTheWall) {

        wall_lat->normal.resize(_dim);
        std::fill(wall_lat->normal.begin(), wall_lat->normal.end(), 0.0);

        int x_pos = wall_lat->x_position;
        int y_pos = wall_lat->y_position;
        int z_pos = wall_lat->z_position;

        std::vector<int> ret = Geometry_.RetrieveLocalInd(x_pos, y_pos, z_pos);

        double x = Geometry_.pointsBox[ret[0]][ret[1]][ret[2]]->x;
        double y = Geometry_.pointsBox[ret[0]][ret[1]][ret[2]]->y;
        double z = Geometry_.pointsBox[ret[0]][ret[1]][ret[2]]->z;

        for (const auto& triangle : Geometry_.localTriangles) {

            // Calculate triangle's normal
            double ux = triangle.x2 - triangle.x1;
            double uy = triangle.y2 - triangle.y1;
            double uz = triangle.z2 - triangle.z1;

            double vx = triangle.x3 - triangle.x1;
            double vy = triangle.y3 - triangle.y1;
            double vz = triangle.z3 - triangle.z1;

            double nx = uy * vz - uz * vy;
            double ny = uz * vx - ux * vz;
            double nz = ux * vy - uy * vx;

            double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
            if (norm > 0) {
                nx /= norm;
                ny /= norm;
                nz /= norm;
            }

            // Calculate area of the triangle (cross product magnitude)
            double area = 0.5 * norm;

            // Compute centroid of the triangle
            double centroidX = (triangle.x1 + triangle.x2 + triangle.x3) / 3.0;
            double centroidY = (triangle.y1 + triangle.y2 + triangle.y3) / 3.0;
            double centroidZ = (triangle.z1 + triangle.z2 + triangle.z3) / 3.0;

            // Compute distance from point to triangle centroid
            double distance = std::sqrt(
                std::pow(x - centroidX, 2) +
                std::pow(y - centroidY, 2) +
                std::pow(z - centroidZ, 2)
            );

            // Weight is the area divided by distance squared (Voronoi area weighting)
            double weight = area / (distance * distance);

            // Accumulate weighted normal
            wall_lat->normal[0] += nx * weight;
            wall_lat->normal[1] += ny * weight;
            wall_lat->normal[2] += nz * weight;
        }

        // Normalize the resulting normal vector
        double norm = std::sqrt(
            wall_lat->normal[0] * wall_lat->normal[0] +
            wall_lat->normal[1] * wall_lat->normal[1] +
            wall_lat->normal[2] * wall_lat->normal[2]
        );

        if (norm > 0) {
            wall_lat->normal[0] /= norm;
            wall_lat->normal[1] /= norm;
            wall_lat->normal[2] /= norm;
        }
    }
}


/*

void SubDomain::calculateBoundaryNormals() {
    for (auto& wall_lat : CertainDepthNode) {

        wall_lat->normal.resize(_dim);
        std::fill(wall_lat->normal.begin(), wall_lat->normal.end(), 0.0);

        int x_pos = wall_lat->x_position;
        int y_pos = wall_lat->y_position;
        int z_pos = wall_lat->z_position;

        std::vector<int> ret = Geometry_.RetrieveLocalInd(x_pos, y_pos, z_pos);

        double x = Geometry_.pointsBox[ret[0]][ret[1]][ret[2]]->x;
        double y = Geometry_.pointsBox[ret[0]][ret[1]][ret[2]]->y;
        double z = Geometry_.pointsBox[ret[0]][ret[1]][ret[2]]->z;

        for (const auto& triangle : Geometry_.localTriangles) {

            // Calculate triangle's normal
            double ux = triangle.x2 - triangle.x1;
            double uy = triangle.y2 - triangle.y1;
            double uz = triangle.z2 - triangle.z1;

            double vx = triangle.x3 - triangle.x1;
            double vy = triangle.y3 - triangle.y1;
            double vz = triangle.z3 - triangle.z1;

            double nx = uy * vz - uz * vy;
            double ny = uz * vx - ux * vz;
            double nz = ux * vy - uy * vx;

            double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
            if (norm > 0) {
                nx /= norm;
                ny /= norm;
                nz /= norm;
            }

            // Calculate angles at the vertices of the triangle
            double angle = CalculateAngle(triangle, x, y, z);

            // Accumulate angle-weighted normal
            wall_lat->normal[0] += nx * angle;
            wall_lat->normal[1] += ny * angle;
            wall_lat->normal[2] += nz * angle;
        }

        // Normalize the resulting normal vector
        double norm = std::sqrt(
            wall_lat->normal[0] * wall_lat->normal[0] +
            wall_lat->normal[1] * wall_lat->normal[1] +
            wall_lat->normal[2] * wall_lat->normal[2]
        );

        if (norm > 0) {
            wall_lat->normal[0] /= norm;
            wall_lat->normal[1] /= norm;
            wall_lat->normal[2] /= norm;
        }
    }
}



double SubDomain::CalculateAngle(const Triangle& triangle, double x, double y, double z) {
    // Determine which vertex of the triangle corresponds to the boundary point (or is closest to it)
    double d1 = std::sqrt(std::pow(triangle.x1 - x, 2) + std::pow(triangle.y1 - y, 2) + std::pow(triangle.z1 - z, 2));
    double d2 = std::sqrt(std::pow(triangle.x2 - x, 2) + std::pow(triangle.y2 - y, 2) + std::pow(triangle.z2 - z, 2));
    double d3 = std::sqrt(std::pow(triangle.x3 - x, 2) + std::pow(triangle.y3 - y, 2) + std::pow(triangle.z3 - z, 2));

    double angle = 0.0;
    if (d1 < d2 && d1 < d3) {
        // The point is closest to vertex 1
        angle = CalculateVertexAngle(triangle.x1, triangle.y1, triangle.z1,
                                     triangle.x2, triangle.y2, triangle.z2,
                                     triangle.x3, triangle.y3, triangle.z3);
    } else if (d2 < d1 && d2 < d3) {
        // The point is closest to vertex 2
        angle = CalculateVertexAngle(triangle.x2, triangle.y2, triangle.z2,
                                     triangle.x1, triangle.y1, triangle.z1,
                                     triangle.x3, triangle.y3, triangle.z3);
    } else {
        // The point is closest to vertex 3
        angle = CalculateVertexAngle(triangle.x3, triangle.y3, triangle.z3,
                                     triangle.x1, triangle.y1, triangle.z1,
                                     triangle.x2, triangle.y2, triangle.z2);
    }

    return angle;
}

// Helper function to calculate the angle at a specific vertex
double SubDomain::CalculateVertexAngle(double vx, double vy, double vz,
                                       double v1x, double v1y, double v1z,
                                       double v2x, double v2y, double v2z) {
    // Vector from vertex to the first other vertex
    double vec1x = v1x - vx;
    double vec1y = v1y - vy;
    double vec1z = v1z - vz;

    // Vector from vertex to the second other vertex
    double vec2x = v2x - vx;
    double vec2y = v2y - vy;
    double vec2z = v2z - vz;

    // Dot product of the two vectors
    double dotProduct = vec1x * vec2x + vec1y * vec2y + vec1z * vec2z;

    // Magnitude of the vectors
    double magnitude1 = std::sqrt(vec1x * vec1x + vec1y * vec1y + vec1z * vec1z);
    double magnitude2 = std::sqrt(vec2x * vec2x + vec2y * vec2y + vec2z * vec2z);

    // Calculate the angle using the dot product and magnitudes
    double angle = std::acos(dotProduct / (magnitude1 * magnitude2));

    return angle;  // Return the angle in radians
}

*/


/*

void SubDomain::calculateBoundaryNormals() {

    
    for (auto& wall_lat : nextToTheWall) {

        wall_lat->normal.resize(_dim);
        double minDistance = std::numeric_limits<double>::infinity();
        Triangle closestTriangle;

        int x_pos = wall_lat->x_position;
        int y_pos = wall_lat->y_position;
        int z_pos = wall_lat->z_position;

        std::vector<int> ret = Geometry_.RetrieveLocalInd(x_pos, y_pos, z_pos);

        double x = Geometry_.pointsBox[ret[0]][ret[1]][ret[2]]->x;
        double y = Geometry_.pointsBox[ret[0]][ret[1]][ret[2]]->y;
        double z = Geometry_.pointsBox[ret[0]][ret[1]][ret[2]]->z;

        for (const auto& triangle : Geometry_.localTriangles) {

            double centroidX = (triangle.x1 + triangle.x2 + triangle.x3) / 3.0;
            double centroidY = (triangle.y1 + triangle.y2 + triangle.y3) / 3.0;
            double centroidZ = (triangle.z1 + triangle.z2 + triangle.z3) / 3.0;

            double distance = std::sqrt(
                std::pow(x - centroidX, 2) +
                std::pow(y - centroidY, 2) +
                std::pow(z - centroidZ, 2)
            );


            if (distance < minDistance) {
                minDistance = distance;
                closestTriangle = triangle;
            }
        }

        wall_lat->normal[0] = closestTriangle.nx;
        wall_lat->normal[1] = closestTriangle.ny;
        wall_lat->normal[2] = closestTriangle.nz;


        //std::cout << "Wall Point (" << x << ", " << y << ", " << z
        //    << ") has normal: (" << wall_lat->normal[0] << ", " << wall_lat->normal[1] << ", " << wall_lat->normal[2] << ")\n";
    }
}

*/


void SubDomain::Initialize()
{
    for (int i = 0; i < my_Nx; i++) {
        for (int j = 0; j < my_Ny; j++) {
            for (int k = 0; k < my_Nz; k++) {

                if (m_lattice[i][j][k] != nullptr) {
                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        m_lattice[i][j][k]->Initialize(m_velSet, dir);
                    }

                }
            }
        }
    }
}

void SubDomain::PrintNonNullNeighbors() {
    for (auto& node : WallBoundaryNode_) {
        if (node != nullptr) {
            int nonNullNeighbors = 0;

            for (int dir = 0; dir < _nLatNodes; dir++) {
                if (node->m_neighbours[dir] != nullptr) {
                    nonNullNeighbors++;
                }
            }

            std::cout << "Node at position (" << node->x_position << ", " << node->y_position << ", " << node->z_position
                << ") has " << nonNullNeighbors << " non-nullptr neighbors." << std::endl;
        }
    }
}

void SubDomain::GetFluidHaloLat() {

    for (int i = 0; i < my_Nx; i++) {
        for (int j = 0; j < my_Ny; j++) {
            for (int k = 0; k < my_Nz; k++) {

                if (m_lattice[i][j][k] != nullptr) {
                    Involved_Lat_.push_back(m_lattice[i][j][k]);

                }
            }
        }
    }
}


void SubDomain::Need_Comm() {

    need_communication.resize(18, 0);

    bool fluidFound[18] = { false }; 

    for (int i = 0; i < my_Nx; i++) {
        for (int j = 0; j < my_Ny; j++) {
            for (int k = 0; k < my_Nz; k++) {

                //left
                if (!fluidFound[0]  && m_lattice[1][j][k]!=nullptr) {
                    need_communication[0] = 1;
                    fluidFound[0] = true;
                }

                //right
                if (!fluidFound[1] && m_lattice[my_Nx - 2][j][k]!=nullptr) {
                    need_communication[1] = 1;
                    fluidFound[1] = true;
                }
                //front
                if (!fluidFound[2] && m_lattice[i][1][k]!=nullptr) {
                    need_communication[2] = 1;
                    fluidFound[2] = true;
                }
                //back
                if (!fluidFound[3] && m_lattice[i][my_Ny - 2][k] !=nullptr) {
                    need_communication[3] = 1;
                    fluidFound[3] = true;
                }
                //up
                if (!fluidFound[4] && m_lattice[i][j][my_Nz - 2]!=nullptr) {
                    need_communication[4] = 1;
                    fluidFound[4] = true;
                }
                //down
                if (!fluidFound[5] && m_lattice[i][j][1]!=nullptr) {
                    need_communication[5] = 1;
                    fluidFound[5] = true;
                }
                //right_up
                if (!fluidFound[6] && m_lattice[my_Nx - 2][j][my_Nz - 2]!=nullptr) {
                    need_communication[6] = 1;
                    fluidFound[6] = true;
                }
                //right_down
                if (!fluidFound[7] && m_lattice[my_Nx - 2][j][1]!=nullptr) {
                    need_communication[7] = 1;
                    fluidFound[7] = true;
                }
                //left_up
                if (!fluidFound[8] && m_lattice[1][j][my_Nz - 2]!=nullptr) {
                    need_communication[8] = 1;
                    fluidFound[8] = true;
                }
                //left_down
                if (!fluidFound[9] && m_lattice[1][j][1]!=nullptr) {
                    need_communication[9] = 1;
                    fluidFound[9] = true;
                }
                //up_front
                if (!fluidFound[10] && m_lattice[i][1][my_Nz - 2]!=nullptr) {
                    need_communication[10] = 1;
                    fluidFound[10] = true;
                }
                //up_back
                if (!fluidFound[11] && m_lattice[i][my_Ny - 2][my_Nz - 2]!=nullptr) {
                    need_communication[11] = 1;
                    fluidFound[11] = true;
                }
                //down_front
                if (!fluidFound[12] && m_lattice[i][1][1]!=nullptr) {
                    need_communication[12] = 1;
                    fluidFound[12] = true;
                }
                //down_back
                if (!fluidFound[13] && m_lattice[i][my_Ny - 2][1]!=nullptr) {
                    need_communication[13] = 1;
                    fluidFound[13] = true;
                }
                //right_front
                if (!fluidFound[14] && m_lattice[my_Nx - 2][1][k]!=nullptr) {
                    need_communication[14] = 1;
                    fluidFound[14] = true;
                }
                //right_back
                if (!fluidFound[15] && m_lattice[my_Nx - 2][my_Ny - 2][k]!=nullptr) {
                    need_communication[15] = 1;
                    fluidFound[15] = true;
                }
                //left_front
                if (!fluidFound[16] && m_lattice[1][1][k]!=nullptr) {
                    need_communication[16] = 1;
                    fluidFound[16] = true;
                }
                //left_back
                if (!fluidFound[17] && m_lattice[1][my_Ny - 2][k]!=nullptr) {
                    need_communication[17] = 1;
                    fluidFound[17] = true;
                }
            }
        }
    }
}


const std::shared_ptr<VelocitySet> SubDomain::GetVelSet()
{
    return m_velSet;
}


void SubDomain::WriteHaloLatticeToVTK() {
    std::string filename = "haloLattice_" + std::to_string(myrank_) + ".vtk";

    std::ofstream vtkFile(filename);

    if (!vtkFile.is_open()) {
        std::cerr << "Error: Unable to open VTK file for writing." << std::endl;
        return;
    }

    int localHaloNodesCount = 0;
    for (int i = 0; i < my_Nx; i++) {
        for (int j = 0; j < my_Ny; j++) {
            for (int k = 0; k < my_Nz; k++) {
                if (m_lattice[i][j][k] != nullptr
                    && m_lattice[i][j][k]->NodeType_ == NodeType::Halo) {
                    localHaloNodesCount++;
                }
            }
        }
    }

    int totalHaloNodesCount = 0;
    MPI_Allreduce(&localHaloNodesCount, &totalHaloNodesCount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    vtkFile << "# vtk DataFile Version 3.0" << std::endl;
    vtkFile << "Halo Lattice Data for Processor " << myrank_ << std::endl;
    vtkFile << "ASCII" << std::endl;
    vtkFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

    std::vector<std::shared_ptr<Node>> haloNodes;

    for (int i = 0; i < my_Nx; i++) {
        for (int j = 0; j < my_Ny; j++) {
            for (int k = 0; k < my_Nz; k++) {
                if (m_lattice[i][j][k] != nullptr
                    && m_lattice[i][j][k]->NodeType_ == NodeType::Halo) {
                    haloNodes.push_back(m_lattice[i][j][k]);
                }
            }
        }
    }

    vtkFile << "POINTS " << haloNodes.size() << " double" << std::endl;
    for (const auto& nodePtr : haloNodes) {
        int i = nodePtr->x_position;
        int j = nodePtr->y_position;
        int k = nodePtr->z_position;

        vtkFile << i << " " << j << " " << k << std::endl;
    }

    vtkFile << "CELLS " << haloNodes.size() << " " << 2 * haloNodes.size() << std::endl;
    for (size_t i = 0; i < haloNodes.size(); ++i) {
        vtkFile << "1 " << i << std::endl;
    }

    vtkFile << "CELL_TYPES " << haloNodes.size() << std::endl;
    for (size_t i = 0; i < haloNodes.size(); ++i) {
        vtkFile << "1" << std::endl;
    }

    vtkFile.close();

    int globalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);
    if (globalRank == 0) {
        std::cout << "Total halo lattice nodes across all processors: " << totalHaloNodesCount << std::endl;
    }
}


void SubDomain::WriteFluidLatticeToVTK() {
    std::string filename = "fluidLattice_" + std::to_string(myrank_) + ".vtk";

    std::ofstream vtkFile(filename);

    if (!vtkFile.is_open()) {
        std::cerr << "Error: Unable to open VTK file for writing." << std::endl;
        return;
    }

    int localFluidNodesCount = 0;
    for (int i = 0; i < my_Nx; i++) {
        for (int j = 0; j < my_Ny; j++) {
            for (int k = 0; k < my_Nz; k++) {
                if (m_lattice[i][j][k] != nullptr
                    && m_lattice[i][j][k]->NodeType_ != NodeType::Halo) {
                    localFluidNodesCount++;
                }
            }
        }
    }

    int totalFluidNodesCount = 0;
    MPI_Allreduce(&localFluidNodesCount, &totalFluidNodesCount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    vtkFile << "# vtk DataFile Version 3.0" << std::endl;
    vtkFile << "Fluid Lattice Data for Processor " << myrank_ << std::endl;
    vtkFile << "ASCII" << std::endl;
    vtkFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

    std::vector<std::shared_ptr<Node>> fluidNodes;

    for (int i = 0; i < my_Nx; i++) {
        for (int j = 0; j < my_Ny; j++) {
            for (int k = 0; k < my_Nz; k++) {
                if (m_lattice[i][j][k] != nullptr
                    && m_lattice[i][j][k]->NodeType_ != NodeType::Halo) {
                    fluidNodes.push_back(m_lattice[i][j][k]);
                }
            }
        }
    }

    vtkFile << "POINTS " << fluidNodes.size() << " double" << std::endl;
    for (const auto& nodePtr : fluidNodes) {

        int i = nodePtr->x_position;
        int j = nodePtr->y_position;
        int k = nodePtr->z_position;

        double x = Geometry_.pointsBox[i][j][k]->x;
        double y = Geometry_.pointsBox[i][j][k]->y;
        double z = Geometry_.pointsBox[i][j][k]->z;

        vtkFile << x << " " << y << " " << z << std::endl;
    }

    vtkFile << "CELLS " << fluidNodes.size() << " " << 2 * fluidNodes.size() << std::endl;
    for (size_t i = 0; i < fluidNodes.size(); ++i) {
        vtkFile << "1 " << i << std::endl;
    }

    vtkFile << "CELL_TYPES " << fluidNodes.size() << std::endl;
    for (size_t i = 0; i < fluidNodes.size(); ++i) {
        vtkFile << "1" << std::endl;
    }

    vtkFile.close();

    int globalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);
    if (globalRank == 0) {
        std::cout << "Total fluid lattice nodes across all processors: " << totalFluidNodesCount << std::endl;
    }
}




