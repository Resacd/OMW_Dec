#include "ComplexGeometry.h"
#include <queue>
#include <unordered_set>
#include <map>
#include <filesystem>


Geometry::Geometry(const std::string& FileName, std::shared_ptr<VelocitySet> velSet, double spacing, char raydirection) {
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    stlFileName = FileName;
    m_velSet = velSet;
    this->spacing = spacing;
    this->raydirection = raydirection;
    Create_cart_comm();
    GetSubBoundingBox();
    Add_Ghost();
    generatePointsInSubBox();
   // writePointsBox(".vtu");
    DistributeTriangles();
    labelPoints();
    modifier();
//    writeLabeledPointsToFile("Labeled.vtu");
    Comm_Ghost();
    OutlierHandler();
    findBoundaryPoints();
//    writePointsToFile(boundaryPoints, "boundary", "boundary");
    classifyBoundaryPoints();
//    writePointsToFile(InletboundaryPoints, "inlet", "inlet");
//    writePointsToFile(OutletboundaryPoints1, "outlet1", "outlet1");
//    writePointsToFile(OutletboundaryPoints2, "outlet2", "outlet2");
//    writePointsToFile(WallboundaryPoints, "wallboundary", "Wallboundary");


}


Geometry::~Geometry() {


}

void Geometry::Create_cart_comm() {
    int const ndims = 3;
    dims[0] = dims[1] = dims[2] = 0;
    int periods[ndims] = { 0, 0, 0 };
    int reorder = 1;
    coords[0] = coords[1] = coords[2] = 0;

    MPI_Dims_create(num_procs, ndims, dims);
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &cart_comm);

    MPI_Cart_shift(cart_comm, 0, 1, &my_left, &my_right);
    MPI_Cart_shift(cart_comm, 1, 1, &my_front, &my_back);
    MPI_Cart_shift(cart_comm, 2, 1, &my_down, &my_up);
    MPI_Comm_rank(cart_comm, &myrank);
    MPI_Cart_coords(cart_comm, myrank, ndims, coords);
    std::vector<int> cart_coords(ndims);
    MPI_Cart_get(cart_comm, ndims, dims, periods, cart_coords.data());

    neighbors.resize(18);

    getNeighborInfo(my_left, neighbors[0]);
    getNeighborInfo(my_right, neighbors[1]);
    getNeighborInfo(my_front, neighbors[2]);
    getNeighborInfo(my_back, neighbors[3]);
    getNeighborInfo(my_up, neighbors[4]);
    getNeighborInfo(my_down, neighbors[5]);

    std::vector<int> right_up_coords = { coords[0] + 1, coords[1], coords[2] + 1 };
    std::vector<int> right_down_coords = { coords[0] + 1, coords[1], coords[2] - 1 };
    std::vector<int> left_up_coords = { coords[0] - 1, coords[1], coords[2] + 1 };
    std::vector<int> left_down_coords = { coords[0] - 1, coords[1], coords[2] - 1 };
    std::vector<int> up_front_coords = { coords[0], coords[1] - 1, coords[2] + 1 };
    std::vector<int> up_back_coords = { coords[0], coords[1] + 1, coords[2] + 1 };
    std::vector<int> down_front_coords = { coords[0], coords[1] - 1, coords[2] - 1 };
    std::vector<int> down_back_coords = { coords[0], coords[1] + 1, coords[2] - 1 };
    std::vector<int> right_front_coords = { coords[0] + 1, coords[1] - 1, coords[2] };
    std::vector<int> right_back_coords = { coords[0] + 1, coords[1] + 1, coords[2] };
    std::vector<int> left_front_coords = { coords[0] - 1, coords[1] - 1, coords[2] };
    std::vector<int> left_back_coords = { coords[0] - 1, coords[1] + 1, coords[2] };

    right_up_ = is_valid_coordinate(right_up_coords) ? get_rank(cart_comm, right_up_coords) : -1;
    right_down_ = is_valid_coordinate(right_down_coords) ? get_rank(cart_comm, right_down_coords) : -1;
    left_up_ = is_valid_coordinate(left_up_coords) ? get_rank(cart_comm, left_up_coords) : -1;
    left_down_ = is_valid_coordinate(left_down_coords) ? get_rank(cart_comm, left_down_coords) : -1;
    up_front_ = is_valid_coordinate(up_front_coords) ? get_rank(cart_comm, up_front_coords) : -1;
    up_back_ = is_valid_coordinate(up_back_coords) ? get_rank(cart_comm, up_back_coords) : -1;
    down_front_ = is_valid_coordinate(down_front_coords) ? get_rank(cart_comm, down_front_coords) : -1;
    down_back_ = is_valid_coordinate(down_back_coords) ? get_rank(cart_comm, down_back_coords) : -1;
    right_front_ = is_valid_coordinate(right_front_coords) ? get_rank(cart_comm, right_front_coords) : -1;
    right_back_ = is_valid_coordinate(right_back_coords) ? get_rank(cart_comm, right_back_coords) : -1;
    left_front_ = is_valid_coordinate(left_front_coords) ? get_rank(cart_comm, left_front_coords) : -1;
    left_back_ = is_valid_coordinate(left_back_coords) ? get_rank(cart_comm, left_back_coords) : -1;

    getNeighborInfo(right_up_, neighbors[6]);
    getNeighborInfo(right_down_, neighbors[7]);
    getNeighborInfo(left_up_, neighbors[8]);
    getNeighborInfo(left_down_, neighbors[9]);
    getNeighborInfo(up_front_, neighbors[10]);
    getNeighborInfo(up_back_, neighbors[11]);
    getNeighborInfo(down_front_, neighbors[12]);
    getNeighborInfo(down_back_, neighbors[13]);
    getNeighborInfo(right_front_, neighbors[14]);
    getNeighborInfo(right_back_, neighbors[15]);
    getNeighborInfo(left_front_, neighbors[16]);
    getNeighborInfo(left_back_, neighbors[17]);
}

void Geometry::getNeighborInfo(int neighbor_rank, Neighbors& neighbor) {
    neighbor.rank = neighbor_rank;
    if (neighbor_rank != MPI_PROC_NULL && neighbor_rank >= 0) {
        MPI_Cart_coords(cart_comm, neighbor_rank, 3, neighbor.coords);
    }
    else {
        neighbor.coords[0] = -1;
        neighbor.coords[1] = -1;
        neighbor.coords[2] = -1;
    }
}

bool Geometry::is_valid_coordinate(std::vector<int>& coords) {
    return coords[0] >= 0 && coords[0] < dims[0] && coords[1] >= 0 && coords[1] < dims[1] && coords[2] >= 0 && coords[2] < dims[2];
}

int Geometry::get_rank(MPI_Comm comm, const std::vector<int>& coords) {
    int rank;
    MPI_Cart_rank(comm, coords.data(), &rank);
    return rank;
}

void Geometry::GetSubBoundingBox() {

//        globalTriangles = Reader.readSTLASCII(stlFileName);
//        localTriangles=globalTriangles;


    if (myrank == 0) {
        globalTriangles = Reader.readSTLASCII(stlFileName);

        Point minPoint = { std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity() };
        Point maxPoint = { -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity() };

        for (const auto& triangle : globalTriangles) {
            minPoint.x = std::min({ minPoint.x, triangle.x1, triangle.x2, triangle.x3 });
            minPoint.y = std::min({ minPoint.y, triangle.y1, triangle.y2 });
            minPoint.z = std::min({ minPoint.z, triangle.z1, triangle.z2 });
            maxPoint.x = std::max({ maxPoint.x, triangle.x1, triangle.x2, triangle.x3 });
            maxPoint.y = std::max({ maxPoint.y, triangle.y1, triangle.y2 });
            maxPoint.z = std::max({ maxPoint.z, triangle.z1, triangle.z2 });
        }

        boundingbox = { minPoint.x, maxPoint.x, minPoint.y, maxPoint.y, minPoint.z, maxPoint.z };
    }

    MPI_Datatype mpiBoundingBoxType;
    MPI_Type_contiguous(6, MPI_DOUBLE, &mpiBoundingBoxType);
    MPI_Type_commit(&mpiBoundingBoxType);

    MPI_Bcast(&boundingbox, 1, mpiBoundingBoxType, 0, MPI_COMM_WORLD);

     nx = std::ceil((boundingbox.max_x - boundingbox.min_x) / spacing) + 1;
     ny = std::ceil((boundingbox.max_y - boundingbox.min_y) / spacing) + 1;
     nz = std::ceil((boundingbox.max_z - boundingbox.min_z) / spacing);

    //if (myrank == 0) {
    //    printf("Processor 0:\n");
    //    printf("  nx = %d\n", nx);
    //    printf("  ny = %d\n", ny);
    //    printf("  nz = %d\n", nz);
    //}

    MPI_Type_free(&mpiBoundingBoxType);


    my_Nx = nx / dims[0];
    my_Ny = ny / dims[1];
    my_Nz = nz / dims[2];

    local_start_x = coords[0] * my_Nx;
    local_start_y = coords[1] * my_Ny;
    local_start_z = coords[2] * my_Nz;
    local_end_x = local_start_x + my_Nx - 1;
    local_end_y = local_start_y + my_Ny - 1;
    local_end_z = local_start_z + my_Nz - 1;

    if (nx % dims[0] > 0 && coords[0] == dims[0] - 1) {
        my_Nx += nx % dims[0];
        local_end_x = local_start_x + my_Nx - 1;
    }
    if (ny % dims[1] > 0 && coords[1] == dims[1] - 1) {
        my_Ny += ny % dims[1];
        local_end_y = local_start_y + my_Ny - 1;
    }
    if (nz % dims[2] > 0 && coords[2] == dims[2] - 1) {
        my_Nz += nz % dims[2];
        local_end_z = local_start_z + my_Nz - 1;
    }

    Point subboxMin = { boundingbox.min_x + local_start_x * spacing, boundingbox.min_y + local_start_y * spacing, boundingbox.min_z + local_start_z * spacing };
    Point subboxMax = { boundingbox.min_x + (local_end_x + 1) * spacing, boundingbox.min_y + (local_end_y + 1) * spacing, boundingbox.min_z + (local_end_z + 1) * spacing };

    subBox = { subboxMin.x, subboxMax.x, subboxMin.y, subboxMax.y, subboxMin.z, subboxMax.z };

    //printf("Process %d:\n"
    //    "  start_x_ = %d\n"
    //    "  end_x_   = %d\n"
    //    "  start_y_ = %d\n"
    //    "  end_y_   = %d\n"
    //    "  start_z_ = %d\n"
    //    "  end_z_   = %d\n",
    //    myrank, local_start_x, local_end_x, local_start_y, local_end_y, local_start_z, local_end_z);
}



void Geometry::Add_Ghost() {

    halo_size = 1;

    if (dims[0] > 1) {
        if (coords[0] == 0) {
            my_Nx += halo_size;
            local_end_x = local_start_x + my_Nx - 1;
        }
        else if (coords[0] == dims[0] - 1) {
            my_Nx += halo_size;
            local_start_x = local_start_x + 2 * halo_size * coords[0] - halo_size;
            local_end_x = local_start_x + my_Nx - 1;
        }
        else {
            my_Nx += 2 * halo_size;
            local_start_x = local_start_x + 2 * halo_size * coords[0] - halo_size;
            local_end_x = local_start_x + my_Nx - 1;
        }
    }

    if (dims[1] > 1) {
        if (coords[1] == 0) {
            my_Ny += halo_size;
            local_end_y = local_start_y + my_Ny - 1;
        }
        else if (coords[1] == dims[1] - 1) {
            my_Ny += halo_size;
            local_start_y = local_start_y + 2 * halo_size * coords[1] - halo_size;
            local_end_y = local_start_y + my_Ny - 1;
        }
        else {
            my_Ny += 2 * halo_size;
            local_start_y = local_start_y + 2 * halo_size * coords[1] - halo_size;
            local_end_y = local_start_y + my_Ny - 1;
        }
    }

    if (dims[2] > 1) {
        if (coords[2] == 0) {
            my_Nz += halo_size;
            local_end_z = local_start_z + my_Nz - 1;
        }
        else if (coords[2] == dims[2] - 1) {
            my_Nz += halo_size;
            local_start_z = local_start_z + 2 * halo_size * coords[2] - halo_size;
            local_end_z = local_start_z + my_Nz - 1;
        }
        else {
            my_Nz += 2 * halo_size;
            local_start_z = local_start_z + 2 * halo_size * coords[2] - halo_size;
            local_end_z = local_start_z + my_Nz - 1;
        }
    }


    //for (int p = 0; p < num_procs; p++) {
    //    if (myrank == p) {
    //        printf("Processor %d:\n", myrank);
    //        printf("my_Nx = %d, my_Ny = %d, my_Nz = %d\n", my_Nx, my_Ny, my_Nz);
    //        printf("start_x = %d, end_x = %d\n", local_start_x, local_end_x);
    //        printf("start_y = %d, end_y = %d\n", local_start_y, local_end_y);
    //        printf("start_z = %d, end_z = %d\n", local_start_z, local_end_z);
    //        fflush(stdout);
    //    }
    //    MPI_Barrier(cart_comm);
    //}

}


bool Geometry::IsGhost(int i, int j, int k) {

    if (my_left != MPI_PROC_NULL && i == 0) {
        return true;
    }

    if (my_right != MPI_PROC_NULL && i == my_Nx - 1) {
        return true;
    }

    if (my_front != MPI_PROC_NULL && j == 0) {
        return true;
    }

    if (my_back != MPI_PROC_NULL && j == my_Ny - 1) {
        return true;
    }

    if (my_up != MPI_PROC_NULL && k == my_Nz - 1) {
        return true;
    }

    if (my_down != MPI_PROC_NULL && k == 0) {
        return true;
    }

    return false;
}

void Geometry::generatePointsInSubBox() {
    pointsBox.resize(my_Nx);
    for (int i = 0; i <my_Nx; ++i) {
        pointsBox[i].resize(my_Ny);
        for (int j = 0; j <my_Ny; ++j) {
            pointsBox[i][j].resize(my_Nz);
        }
    }

    for (int k = 0; k <my_Nz; ++k) {
        for (int j = 0; j < my_Ny; ++j) {
            for (int i = 0; i < my_Nx; ++i) {
                auto point = std::make_shared<Point>();

                if (IsGhost(i, j, k)) {
                    point->PointBoxType = PointBoxType::ghost;
                }
                else {

                    std::vector<int> l2g = Local2Global(i, j, k);
                    double x = boundingbox.min_x + l2g[0] * spacing;
                    double y = boundingbox.min_y + l2g[1] * spacing;
                    double z = boundingbox.min_z + l2g[2] * spacing;
                    point->x = x;
                    point->y = y;
                    point->z = z;
                    point->PointBoxType = PointBoxType::real;
                }

                pointsBox[i][j][k] = point;
            }
        }
    }
}

bool Geometry::IsInBox(int coord_x, int coord_y, int coord_z) {

    if (coord_x >= 0 && coord_x < pointsBox.size()
        && coord_y >= 0 && coord_y < pointsBox[0].size()
        && coord_z >= 0 && coord_z < pointsBox[0][0].size())
        return true;
    else
        return false;

}


std::vector<int>Geometry::Local2Global(int i, int j, int k) {

    std::vector<int> L2G;
    int global_coordX = i;
    int global_coordY = j;
    int global_coordZ = k;

    if (dims[0] > 1) {
        if (coords[0] == 0) {
            global_coordX = i;
        }
        else if (coords[0] == dims[0] - 1) {

            global_coordX = local_start_x + i - 2 * coords[0] * halo_size;
        }
        else {

            global_coordX = local_start_x + i - 2 * coords[0] * halo_size;

        }
    }


    if (dims[1] > 1) {
        if (coords[1] == 0) {
            ;
            global_coordY = j;
        }
        else if (coords[1] == dims[1] - 1) {

            global_coordY = local_start_y + j - 2 * coords[1] * halo_size;
        }
        else {

            global_coordY = local_start_y + j - 2 * coords[1] * halo_size;

        }
    }

    if (dims[2] > 1) {
        if (coords[2] == 0) {
            ;
            global_coordZ = k;
        }
        else if (coords[2] == dims[2] - 1) {

            global_coordZ = local_start_z + k - 2 * coords[2] * halo_size;
        }
        else {

            global_coordZ = local_start_z + k - 2 * coords[2] * halo_size;

        }
    }

    L2G = { global_coordX,global_coordY,global_coordZ };
    return L2G;

}


std::vector<int>Geometry::RetrieveLocalInd(int i, int j, int k) {

    std::vector<int> Retrieved;
    int local_coordX = i;
    int local_coordY = j;
    int local_coordZ = k;

    if (dims[0] > 1) {
        if (coords[0] == 0) {
            local_coordX = i;
        }
        else if (coords[0] == dims[0] - 1) {

            local_coordX = i - halo_size;
        }
        else {

            local_coordX = i - halo_size;

        }
    }


    if (dims[1] > 1) {
        if (coords[1] == 0) {
            ;
            local_coordY = j;
        }
        else if (coords[1] == dims[1] - 1) {

            local_coordY = j - halo_size;
        }
        else {

            local_coordY = j - halo_size;

        }
    }

    if (dims[2] > 1) {
        if (coords[2] == 0) {
            ;
            local_coordZ = k;
        }
        else if (coords[2] == dims[2] - 1) {

            local_coordZ = k - halo_size;
        }
        else {

            local_coordZ = k - halo_size;

        }
    }

    Retrieved = { local_coordX,local_coordY,local_coordZ };
    return Retrieved;

}



bool pointInBox(double x, double y, double z, const BoundingBox& box, char raydirection) {
    switch (raydirection) {
        case 'x':
            return (y >= box.min_y && y <= box.max_y) && (z >= box.min_z && z <= box.max_z);
        case 'y':
            return (x >= box.min_x && x <= box.max_x) && (z >= box.min_z && z <= box.max_z);
        case 'z':
            return (x >= box.min_x && x <= box.max_x) && (y >= box.min_y && y <= box.max_y);
        default:
            return false;
    }
}

bool edgeIntersectsBox(double x1, double y1, double z1, double x2, double y2, double z2, const BoundingBox& box) {
    // This function should check if the edge between (x1, y1, z1) and (x2, y2, z2) intersects the bounding box
    // This is a simplified version, a more comprehensive intersection check would be needed
    return (std::min(x1, x2) <= box.max_x && std::max(x1, x2) >= box.min_x) &&
           (std::min(y1, y2) <= box.max_y && std::max(y1, y2) >= box.min_y) &&
           (std::min(z1, z2) <= box.max_z && std::max(z1, z2) >= box.min_z);
}

bool planeIntersectsTriangle(double plane_d, double plane_nx, double plane_ny, double plane_nz,
                             const Triangle& tri) {
    // This function should check if the plane intersects the triangle
    // A comprehensive intersection test would be required here
    double dist1 = plane_nx * tri.x1 + plane_ny * tri.y1 + plane_nz * tri.z1 - plane_d;
    double dist2 = plane_nx * tri.x2 + plane_ny * tri.y2 + plane_nz * tri.z2 - plane_d;
    double dist3 = plane_nx * tri.x3 + plane_ny * tri.y3 + plane_nz * tri.z3 - plane_d;
    
    // Check if the triangle is straddling the plane
    return (dist1 * dist2 <= 0) || (dist2 * dist3 <= 0) || (dist1 * dist3 <= 0);
}

bool Geometry::DoesTriangleIntersectBox(const Triangle& tri, const BoundingBox& box) {
    // Check if any of the triangle's vertices are inside the box
    if (pointInBox(tri.x1, tri.y1, tri.z1, box, this->raydirection) ||
        pointInBox(tri.x2, tri.y2, tri.z2, box, this->raydirection) ||
        pointInBox(tri.x3, tri.y3, tri.z3, box, this->raydirection)) {
        return true;
    }

    // Check if any of the triangle's edges intersect the box
    if (edgeIntersectsBox(tri.x1, tri.y1, tri.z1, tri.x2, tri.y2, tri.z2, box) ||
        edgeIntersectsBox(tri.x2, tri.y2, tri.z2, tri.x3, tri.y3, tri.z3, box) ||
        edgeIntersectsBox(tri.x3, tri.y3, tri.z3, tri.x1, tri.y1, tri.z1, box)) {
        return true;
    }

    // Check if the triangle intersects any of the box's planes
    if (planeIntersectsTriangle(box.min_x, 1, 0, 0, tri) ||
        planeIntersectsTriangle(box.max_x, -1, 0, 0, tri) ||
        planeIntersectsTriangle(box.min_y, 0, 1, 0, tri) ||
        planeIntersectsTriangle(box.max_y, 0, -1, 0, tri) ||
        planeIntersectsTriangle(box.min_z, 0, 0, 1, tri) ||
        planeIntersectsTriangle(box.max_z, 0, 0, -1, tri)) {
        return true;
    }

    return false;
}




void Geometry::DistributeTriangles() {

    std::vector<BoundingBox> allSubBoxes(num_procs);
    MPI_Allgather(&subBox, sizeof(BoundingBox), MPI_BYTE, allSubBoxes.data(), sizeof(BoundingBox), MPI_BYTE, MPI_COMM_WORLD);

    std::vector<std::vector<Triangle>> sendBuffers(num_procs);
    if (myrank == 0) {
        for (const auto& triangle : globalTriangles) {
            for (int i = 0; i < num_procs; ++i) {
                if (DoesTriangleIntersectBox(triangle,allSubBoxes[i])) {
                    sendBuffers[i].push_back(triangle);
                }
            }
        }

        localTriangles = sendBuffers[0];
        for (int i = 1; i < num_procs; ++i) {
            int numTriangles = sendBuffers[i].size();
            MPI_Send(&numTriangles, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(sendBuffers[i].data(), numTriangles * sizeof(Triangle), MPI_BYTE, i, 0, MPI_COMM_WORLD);
        }
    }
    else {
        int numTriangles;
        MPI_Recv(&numTriangles, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        localTriangles.resize(numTriangles);
        MPI_Recv(localTriangles.data(), numTriangles * sizeof(Triangle), MPI_BYTE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

       // std::cout << "Processor " << myrank << " received " << localTriangles.size() << " triangles." << std::endl;
    }

    if (myrank == 0) {
       // std::cout << "Processor " << myrank << " has " << sendBuffers[0].size() << " triangles." << std::endl;
    }
}


bool Geometry::doLinesIntersect(const Line& line1, const Line& line2) {

    double x1 = line1.start.x;
    double y1 = line1.start.y;
    double z1 = line1.start.z;
    double x2 = line1.end.x;
    double y2 = line1.end.y;
    double z2 = line1.end.z;

    double x3 = line2.start.x;
    double y3 = line2.start.y;
    double z3 = line2.start.z;
    double x4 = line2.end.x;
    double y4 = line2.end.y;
    double z4 = line2.end.z;

    double ua = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / ((y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1));
    double ub = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / ((y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1));

    const double epsilon = 1e-7;
    if (ua >= -epsilon && ua <= 1 + epsilon && ub >= -epsilon && ub <= 1 + epsilon) {
        return true;
    }

    return false;
}


bool Geometry::isPointInsideTriangle(double x, double y, const Triangle& triangle) {

    double x1 = triangle.x1;
    double y1 = triangle.y1;
    double x2 = triangle.x2;
    double y2 = triangle.y2;
    double x3 = triangle.x3;
    double y3 = triangle.y3;

    double d1 = ((x - x2) * (y1 - y2) - (y - y2) * (x1 - x2));
    double d2 = ((x - x3) * (y2 - y3) - (y - y3) * (x2 - x3));
    double d3 = ((x - x1) * (y3 - y1) - (y - y1) * (x3 - x1));

    bool hasNegative = (d1 < 0) || (d2 < 0) || (d3 < 0);
    bool hasPositive = (d1 > 0) || (d2 > 0) || (d3 > 0);

    return !(hasNegative && hasPositive);
}


int Geometry::countIntersectinsTriangle(Line& line, std::vector<Triangle>& posintersectedtri) {
    const double threshold = 1e-7;
    int count = 0;

    for (const auto& triangle : posintersectedtri) {

        Line edge1, edge2, edge3;
        edge1.start = { triangle.x1, triangle.y1, triangle.z1 };
        edge1.end = { triangle.x2, triangle.y2, triangle.z2 };
        edge2.start = { triangle.x2, triangle.y2, triangle.z2 };
        edge2.end = { triangle.x3, triangle.y3, triangle.z3 };
        edge3.start = { triangle.x3, triangle.y3, triangle.z3 };
        edge3.end = { triangle.x1, triangle.y1, triangle.z1 };

        //if (doLinesIntersect(line, edge1) || doLinesIntersect(line, edge2) || doLinesIntersect(line, edge3)) {
        //    count++;
        //    continue;
        //}

        double px = line.start.x;
        double py = line.start.y;
        double pz = line.start.z;

        double nx = triangle.nx;
        double ny = triangle.ny;
        double nz = triangle.nz;

        double qx = triangle.x1;
        double qy = triangle.y1;
        double qz = triangle.z1;

        double t = ((qx - px) * nx + (qy - py) * ny + (qz - pz) * nz) /
            ((line.end.x - line.start.x) * nx + (line.end.y - line.start.y) * ny + (line.end.z - line.start.z) * nz);

        if (t >= -threshold && t <= 1 + threshold) {
            double ix = px + t * (line.end.x - line.start.x);
            double iy = py + t * (line.end.y - line.start.y);
            double iz = pz + t * (line.end.z - line.start.z);

            if (isPointInsideTriangle(ix, iy, triangle) && t >= -threshold && t <= 1 + threshold) {
                count++;
            }
        }
    }

    return count;
}


void Geometry::adjustLineCoordinates(Line& line, BoundingBox& boundingbox, double rayOffset, int directionSwitch) {
    switch (directionSwitch) {
    case 0:  // positive z-direction
        if (line.start.x == boundingbox.min_x) {
            line.start.x += rayOffset;
            line.end.x += rayOffset;
        }

        if (line.start.x == boundingbox.max_x) {
            line.start.x -= rayOffset;
            line.end.x -= rayOffset;
        }

        if (line.start.y == boundingbox.min_y) {
            line.start.y += rayOffset;
            line.end.y += rayOffset;
        }

        if (line.start.y == boundingbox.max_y) {
            line.start.y -= rayOffset;
            line.end.y -= rayOffset;
        }

        if (line.start.z == boundingbox.min_z) {
            line.start.z += rayOffset;
        }

        if (line.end.z == boundingbox.max_z) {
            line.end.z -= rayOffset;
        }

        break;

    case 1:  // negative z-direction
        if (line.start.x == boundingbox.min_x) {
            line.start.x += rayOffset;
            line.end.x += rayOffset;
        }

        if (line.start.x == boundingbox.max_x) {
            line.start.x -= rayOffset;
            line.end.x -= rayOffset;
        }

        if (line.start.y == boundingbox.min_y) {
            line.start.y += rayOffset;
            line.end.y += rayOffset;
        }

        if (line.start.y == boundingbox.max_y) {
            line.start.y -= rayOffset;
            line.end.y -= rayOffset;
        }


        if (line.start.z == boundingbox.min_z) {
            line.start.z -= rayOffset;
        }

        break;

    case 2:  // positive x-direction

        if (line.start.z == boundingbox.min_z) {
            line.start.z += rayOffset;
            line.end.z += rayOffset;
        }

        if (line.start.z == boundingbox.max_z) {
            line.start.z -= rayOffset;
            line.end.z -= rayOffset;
        }

        if (line.start.y == boundingbox.min_y) {
            line.start.y += rayOffset;
            line.end.y += rayOffset;
        }

        if (line.start.y == boundingbox.max_y) {
            line.start.y -= rayOffset;
            line.end.y -= rayOffset;
        }

        if (line.start.x == boundingbox.min_x) {
            line.start.x += rayOffset;
            line.end.x += rayOffset;
        }

        if (line.start.x == boundingbox.max_x) {
            line.start.x -= rayOffset;
            line.end.x -= rayOffset;
        }

        break;


    case 3:  // negative x-direction


    case 4:  // positive y-direction
    case 5:  // negative y-direction
        if (line.start.x == boundingbox.min_x) {
            line.start.x += rayOffset;
            line.end.x += rayOffset;
        }

        if (line.start.x == boundingbox.max_x) {
            line.start.x -= rayOffset;
            line.end.x -= rayOffset;
        }

        if (line.start.z == boundingbox.min_z) {
            line.start.z += rayOffset;
            line.end.z += rayOffset;
        }

        if (line.start.z == boundingbox.max_z) {
            line.start.z -= rayOffset;
            line.end.z -= rayOffset;
        }
        break;
    }
}


void Geometry::labelPoints() {
    const double rayOffset = 0.02;
    const double tolerance = 1e-6;
    int directionSwitch;
    int local_intersected_triangles;

    for (int i = 0; i < pointsBox.size(); i++) {
        for (int j = 0; j < pointsBox[i].size(); j++) {
            for (int k = 0; k < pointsBox[i][j].size(); k++) {

                if (pointsBox[i][j][k]->PointBoxType == PointBoxType::real) {

                    Line line;
                    line.start.x = pointsBox[i][j][k]->x;
                    line.start.y = pointsBox[i][j][k]->y;
                    line.start.z = pointsBox[i][j][k]->z;

                    directionSwitch = 0;
                    local_intersected_triangles = 0;
                    bool foundNonZeroIntersections = false;

                    switch (directionSwitch) {
                    case 0:  // positive z-direction
                        rayLength = boundingbox.max_z - boundingbox.min_z;
                        line.end.x = pointsBox[i][j][k]->x;
                        line.end.y = pointsBox[i][j][k]->y;
                        line.end.z = pointsBox[i][j][k]->z + rayLength + 0.1;
                        break;
                    case 1:  // negative z-direction
                        rayLength = boundingbox.max_z - boundingbox.min_z;
                        line.end.x = pointsBox[i][j][k]->x;
                        line.end.y = pointsBox[i][j][k]->y;
                        line.end.z = pointsBox[i][j][k]->z - (rayLength + 0.1);
                        break;
                    case 2:  // positive x-direction
                        rayLength = boundingbox.max_x - boundingbox.min_x;
                        line.end.x = pointsBox[i][j][k]->x + rayLength + 0.1;
                        line.end.y = pointsBox[i][j][k]->y;
                        line.end.z = pointsBox[i][j][k]->z;
                        break;
                    case 3:  // negative x-direction
                        rayLength = boundingbox.max_x - boundingbox.min_x;
                        line.end.x = pointsBox[i][j][k]->x - (rayLength + 0.1);
                        line.end.y = pointsBox[i][j][k]->y;
                        line.end.z = pointsBox[i][j][k]->z;
                        break;
                    case 4:  // positive y-direction
                        rayLength = boundingbox.max_y - boundingbox.min_y;
                        line.end.x = pointsBox[i][j][k]->x;
                        line.end.y = pointsBox[i][j][k]->y + rayLength + 0.1;
                        line.end.z = pointsBox[i][j][k]->z;
                        break;
                    case 5:  // negative y-direction
                        rayLength = boundingbox.max_y - boundingbox.min_y;
                        line.end.x = pointsBox[i][j][k]->x;
                        line.end.y = pointsBox[i][j][k]->y - (rayLength + 0.1);
                        line.end.z = pointsBox[i][j][k]->z;
                        break;
                    }

                    adjustLineCoordinates(line, boundingbox, rayOffset, directionSwitch);

                    local_intersected_triangles = countIntersectinsTriangle(line, localTriangles);


                    if (local_intersected_triangles % 2 == 1) {
                        pointsBox[i][j][k]->NodeType = NodeType::Fluid;
                    }
                    else {
                        pointsBox[i][j][k]->NodeType = NodeType::NotFluid;
                    }
                }
            }
        }
    }
}


/*
void Geometry::OutlierHandler() {
    std::vector<std::tuple<int, int, int>> pointsToRelabel;

    for (int i = 0; i < pointsBox.size(); i++) {
        for (int j = 0; j < pointsBox[i].size(); j++) {
            for (int k = 0; k < pointsBox[i][j].size(); k++) {
                
                
                if (pointsBox[i][j][k]->NodeType == NodeType::Fluid && pointsBox[i][j][k]->PointBoxType == PointBoxType::real) {
                    bool isSingular = false;

                    
                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        int opp_dir = m_velSet->OppositeDirection(dir);

                        int x_dir = i + m_velSet->GetDirection(dir)[0];
                        int y_dir = j + m_velSet->GetDirection(dir)[1];
                        int z_dir = k + m_velSet->GetDirection(dir)[2];

                        int x_opp_dir = i + m_velSet->GetDirection(opp_dir)[0];
                        int y_opp_dir = j + m_velSet->GetDirection(opp_dir)[1];
                        int z_opp_dir = k + m_velSet->GetDirection(opp_dir)[2];

                        bool dirExists = IsInBox(x_dir, y_dir, z_dir);
                        bool oppDirExists = IsInBox(x_opp_dir, y_opp_dir, z_opp_dir);

                        bool dirIsNotFluid = dirExists && pointsBox[x_dir][y_dir][z_dir]->NodeType == NodeType::NotFluid;
                        bool oppDirIsNotFluid = oppDirExists && pointsBox[x_opp_dir][y_opp_dir][z_opp_dir]->NodeType == NodeType::NotFluid;

                        
                        if ((!dirExists || dirIsNotFluid) && (!oppDirExists || oppDirIsNotFluid)) {
                            isSingular = true;
                            break;
                        }
                    }

                    if (isSingular) {
                        pointsToRelabel.push_back({i, j, k});  
                    }
                }
            }
        }
    }


    for (const auto& point : pointsToRelabel) {
        int i = std::get<0>(point);
        int j = std::get<1>(point);
        int k = std::get<2>(point);

        pointsBox[i][j][k]->NodeType = NodeType::NotFluid;
    }
}


*/




void Geometry::OutlierHandler() {
    std::vector<std::tuple<int, int, int>> pointsToRelabel;

    for (int i = 0; i < pointsBox.size(); i++) {
        for (int j = 0; j < pointsBox[i].size(); j++) {
            for (int k = 0; k < pointsBox[i][j].size(); k++) {
                
                if (pointsBox[i][j][k]->NodeType == NodeType::Fluid && pointsBox[i][j][k]->PointBoxType == PointBoxType::real) {
                    bool isSingular = false;
                    int fluidNeighborCount = 0;

                      for (int dir = 0; dir < 18; dir++) {
                        int x_dir = i + m_velSet->GetDirection(dir)[0];
                        int y_dir = j + m_velSet->GetDirection(dir)[1];
                        int z_dir = k + m_velSet->GetDirection(dir)[2];

                        int opp_dir = m_velSet->OppositeDirection(dir);
                        int x_opp_dir = i + m_velSet->GetDirection(opp_dir)[0];
                        int y_opp_dir = j + m_velSet->GetDirection(opp_dir)[1];
                        int z_opp_dir = k + m_velSet->GetDirection(opp_dir)[2];

                        bool dirExists = IsInBox(x_dir, y_dir, z_dir);
                        bool oppDirExists = IsInBox(x_opp_dir, y_opp_dir, z_opp_dir);

                        bool dirIsNotFluid = dirExists && pointsBox[x_dir][y_dir][z_dir]->NodeType == NodeType::NotFluid;
                        bool oppDirIsNotFluid = oppDirExists && pointsBox[x_opp_dir][y_opp_dir][z_opp_dir]->NodeType == NodeType::NotFluid;

                        if ((!dirExists || dirIsNotFluid) && (!oppDirExists || oppDirIsNotFluid)) {
                            isSingular = true;
                            break;
                        }


                        if (dirExists && pointsBox[x_dir][y_dir][z_dir]->NodeType == NodeType::Fluid) {
                            fluidNeighborCount++;
                        }
                    }


                    if (fluidNeighborCount < 6) {
                        isSingular = true;
                    }

                    if (isSingular) {
                        pointsToRelabel.push_back({i, j, k});
                    }
                }
            }
        }
    }


    for (const auto& point : pointsToRelabel) {
        int i = std::get<0>(point);
        int j = std::get<1>(point);
        int k = std::get<2>(point);

        pointsBox[i][j][k]->NodeType = NodeType::NotFluid;
    }
}




void Geometry::modifier() {
    const double tolerance = 1e-2;  

    for (int i = 0; i < pointsBox.size(); i++) {
        for (int j = 0; j < pointsBox[i].size(); j++) {
            for (int k = 0; k < pointsBox[i][j].size(); k++) {

                if (std::abs(pointsBox[i][j][k]->x - (-20.6236)) < tolerance &&
                    std::abs(pointsBox[i][j][k]->y - (-229.439)) < tolerance &&
                    pointsBox[i][j][k]->z <= -121.565) {

                    pointsBox[i][j][k]->NodeType = NodeType::NotFluid;
                }
            }
        }
    }
}


void Geometry::Comm_Ghost() {
    int tag_xp = 0;
    int tag_xm = 1;
    int tag_yp = 2;
    int tag_ym = 3;
    int tag_zp = 4;
    int tag_zm = 5;

    int mpi_error;
    MPI_Status status;

    auto pack_point = [](std::shared_ptr<Point> point) {
        return std::vector<double>{ point->x, point->y, point->z, static_cast<double>(point->NodeType) };
    };

    auto unpack_point = [](const std::vector<double>& data, size_t& index) {
        double x = data[index++];
        double y = data[index++];
        double z = data[index++];
        NodeType nodeType = static_cast<NodeType>(data[index++]);
        return std::make_tuple(x, y, z, nodeType);
    };

    if (my_left != MPI_PROC_NULL) {
        std::vector<double> send_data;
        for (int k = 0; k < my_Nz; k++) {
            for (int j = 0; j < my_Ny; j++) {
                auto point_data = pack_point(pointsBox[1][j][k]);
                send_data.insert(send_data.end(), point_data.begin(), point_data.end());
            }
        }

        std::vector<double> recv_data(send_data.size());
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, my_left, tag_xm,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, my_left, tag_xp,
            MPI_COMM_WORLD, &status);

        size_t index = 0;
        for (int k = 0; k < my_Nz; k++) {
            for (int j = 0; j < my_Ny; j++) {
                double x, y, z;
                NodeType nodeType;
                std::tie(x, y, z, nodeType) = unpack_point(recv_data, index);
                pointsBox[0][j][k]->x = x;
                pointsBox[0][j][k]->y = y;
                pointsBox[0][j][k]->z = z;
                pointsBox[0][j][k]->NodeType = nodeType;
            }
        }
    }

    // my_right
    if (my_right != MPI_PROC_NULL) {
        std::vector<double> send_data;
        for (int k = 0; k < my_Nz; k++) {
            for (int j = 0; j < my_Ny; j++) {
                auto point_data = pack_point(pointsBox[my_Nx - 2][j][k]);
                send_data.insert(send_data.end(), point_data.begin(), point_data.end());
            }
        }

        std::vector<double> recv_data(send_data.size());
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, my_right, tag_xp,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, my_right, tag_xm,
            MPI_COMM_WORLD, &status);

        size_t index = 0;
        for (int k = 0; k < my_Nz; k++) {
            for (int j = 0; j < my_Ny; j++) {
                double x, y, z;
                NodeType nodeType;
                std::tie(x, y, z, nodeType) = unpack_point(recv_data, index);
                pointsBox[my_Nx - 1][j][k]->x = x;
                pointsBox[my_Nx - 1][j][k]->y = y;
                pointsBox[my_Nx - 1][j][k]->z = z;
                pointsBox[my_Nx - 1][j][k]->NodeType = nodeType;
            }
        }
    }

    // my_front
    if (my_front != MPI_PROC_NULL) {
        std::vector<double> send_data;
        for (int k = 0; k < my_Nz; k++) {
            for (int i = 0; i < my_Nx; i++) {
                auto point_data = pack_point(pointsBox[i][1][k]);
                send_data.insert(send_data.end(), point_data.begin(), point_data.end());
            }
        }

        std::vector<double> recv_data(send_data.size());
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, my_front, tag_ym,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, my_front, tag_yp,
            MPI_COMM_WORLD, &status);

        size_t index = 0;
        for (int k = 0; k < my_Nz; k++) {
            for (int i = 0; i < my_Nx; i++) {
                double x, y, z;
                NodeType nodeType;
                std::tie(x, y, z, nodeType) = unpack_point(recv_data, index);
                pointsBox[i][0][k]->x = x;
                pointsBox[i][0][k]->y = y;
                pointsBox[i][0][k]->z = z;
                pointsBox[i][0][k]->NodeType = nodeType;
            }
        }
    }

    // my_back
    if (my_back != MPI_PROC_NULL) {
        std::vector<double> send_data;
        for (int k = 0; k < my_Nz; k++) {
            for (int i = 0; i < my_Nx; i++) {
                auto point_data = pack_point(pointsBox[i][my_Ny - 2][k]);
                send_data.insert(send_data.end(), point_data.begin(), point_data.end());
            }
        }

        std::vector<double> recv_data(send_data.size());
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, my_back, tag_yp,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, my_back, tag_ym,
            MPI_COMM_WORLD, &status);

        size_t index = 0;
        for (int k = 0; k < my_Nz; k++) {
            for (int i = 0; i < my_Nx; i++) {
                double x, y, z;
                NodeType nodeType;
                std::tie(x, y, z, nodeType) = unpack_point(recv_data, index);
                pointsBox[i][my_Ny - 1][k]->x = x;
                pointsBox[i][my_Ny - 1][k]->y = y;
                pointsBox[i][my_Ny - 1][k]->z = z;
                pointsBox[i][my_Ny - 1][k]->NodeType = nodeType;
            }
        }
    }

    // my_down
    if (my_down != MPI_PROC_NULL) {
        std::vector<double> send_data;
        for (int j = 0; j < my_Ny; j++) {
            for (int i = 0; i < my_Nx; i++) {
                auto point_data = pack_point(pointsBox[i][j][1]);
                send_data.insert(send_data.end(), point_data.begin(), point_data.end());
            }
        }

        std::vector<double> recv_data(send_data.size());
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, my_down, tag_zm,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, my_down, tag_zp,
            MPI_COMM_WORLD, &status);

        size_t index = 0;
        for (int j = 0; j < my_Ny; j++) {
            for (int i = 0; i < my_Nx; i++) {
                double x, y, z;
                NodeType nodeType;
                std::tie(x, y, z, nodeType) = unpack_point(recv_data, index);
                pointsBox[i][j][0]->x = x;
                pointsBox[i][j][0]->y = y;
                pointsBox[i][j][0]->z = z;
                pointsBox[i][j][0]->NodeType = nodeType;
            }
        }
    }

    // my_up
    if (my_up != MPI_PROC_NULL) {
        std::vector<double> send_data;
        for (int j = 0; j < my_Ny; j++) {
            for (int i = 0; i < my_Nx; i++) {
                auto point_data = pack_point(pointsBox[i][j][my_Nz - 2]);
                send_data.insert(send_data.end(), point_data.begin(), point_data.end());
            }
        }

        std::vector<double> recv_data(send_data.size());
        mpi_error = MPI_Sendrecv(send_data.data(), send_data.size(), MPI_DOUBLE, my_up, tag_zp,
            recv_data.data(), recv_data.size(), MPI_DOUBLE, my_up, tag_zm,
            MPI_COMM_WORLD, &status);

        size_t index = 0;
        for (int j = 0; j < my_Ny; j++) {
            for (int i = 0; i < my_Nx; i++) {
                double x, y, z;
                NodeType nodeType;
                std::tie(x, y, z, nodeType) = unpack_point(recv_data, index);
                pointsBox[i][j][my_Nz - 1]->x = x;
                pointsBox[i][j][my_Nz - 1]->y = y;
                pointsBox[i][j][my_Nz - 1]->z = z;
                pointsBox[i][j][my_Nz - 1]->NodeType = nodeType;
            }
        }
    }
}



std::vector<Point> Geometry::findBoundaryPoints() {

    for (int i = 0; i < pointsBox.size(); i++) {
        for (int j = 0; j < pointsBox[i].size(); j++) {
            for (int k = 0; k < pointsBox[i][j].size(); k++) {
                std::shared_ptr<Point> point = pointsBox[i][j][k];

                if (point->NodeType==NodeType::Fluid && point->PointBoxType!=PointBoxType::ghost) {
                    bool isBoundaryPoint = false;

                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        int x_neighbor = i + m_velSet->GetDirection(dir)[0];
                        int y_neighbor = j + m_velSet->GetDirection(dir)[1];
                        int z_neighbor = k + m_velSet->GetDirection(dir)[2];

                        if (x_neighbor >= 0 && x_neighbor < pointsBox.size() &&
                            y_neighbor >= 0 && y_neighbor < pointsBox[i].size() &&
                            z_neighbor >= 0 && z_neighbor < pointsBox[i][j].size()) {
                            if (pointsBox[x_neighbor][y_neighbor][z_neighbor]->NodeType==NodeType::NotFluid) {
                                isBoundaryPoint = true;
                                break;  
                            }
                        }
                        else {
                            isBoundaryPoint = true;
                            break;
                        }
                    }

                    if (isBoundaryPoint) {
                        boundaryPoints.push_back(*point);
                    }
                }
            }
        }
    }
//    std::cout << "Processor " << myrank << " has " << boundaryPoints.size() << " boundary points." << std::endl;


    return boundaryPoints;
}

void Geometry::classifyBoundaryPoints() {
    long double tolerance = 1e-2;
    std::unordered_set<int> inOutIndices;

    for (int i = 0; i < boundaryPoints.size(); i++) {
        long double x = boundaryPoints[i].x;
        long double y = boundaryPoints[i].y;
        long double z = boundaryPoints[i].z;

        long double zb_inlet = -144.965L;
        if (std::abs(z - zb_inlet) <= tolerance) {
            InletboundaryPoints.push_back(boundaryPoints[i]);
            inOutIndices.insert(i);
        }

        long double xb_outlet1 = -20.7736L;
        if (std::abs(x - xb_outlet1) <= tolerance) {
            OutletboundaryPoints1.push_back(boundaryPoints[i]);
            inOutIndices.insert(i);
        }

        long double zb_outlet2 = -116.965L;
        if ((std::abs(z - zb_outlet2) <= tolerance) &&
            (3.9 < x && x < 6.7) && (-233.0 < y && y < -230.0)) {
            OutletboundaryPoints2.push_back(boundaryPoints[i]);
            inOutIndices.insert(i);
        }
    }

    std::vector<Point> wallPoints;
    for (int i = 0; i < boundaryPoints.size(); ++i) {
        if (inOutIndices.find(i) == inOutIndices.end()) {
            wallPoints.push_back(boundaryPoints[i]);
        }
    }

    this->WallboundaryPoints = wallPoints;
}


void Geometry::writePointsBox(const std::string& filename) {


#ifdef _WIN32
    std::filesystem::create_directory("pointsBox");
#else 
    mkdir("pointsBox", 0777);
#endif

    std::string filepath = "pointsBox/points_" + std::to_string(myrank) + "_" + filename + ".vtk";

    std::ofstream vtkFile(filepath);

    if (!vtkFile.is_open()) {
        std::cerr << "Error: Unable to open VTK file for writing." << std::endl;
        return;
    }

    vtkFile << "# vtk DataFile Version 3.0" << std::endl;
    vtkFile << "Fluid Points Data" << std::endl;
    vtkFile << "ASCII" << std::endl;
    vtkFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

    std::vector<Point> fluidPoints;
    std::vector<int> indicesI, indicesJ, indicesK;

    for (size_t i = 0; i < pointsBox.size(); ++i) {
        for (size_t j = 0; j < pointsBox[i].size(); ++j) {
            for (size_t k = 0; k < pointsBox[i][j].size(); ++k) {
                const auto& point = pointsBox[i][j][k];
                if (point && point->PointBoxType == PointBoxType::real) {
                    fluidPoints.push_back(*point);

                    std::vector<int> l2g = Local2Global(i, j, k);
                    indicesI.push_back(static_cast<int>(l2g[0]));
                    indicesJ.push_back(static_cast<int>(l2g[1]));
                    indicesK.push_back(static_cast<int>(l2g[2]));
                }
            }
        }
    }

    vtkFile << "POINTS " << fluidPoints.size() << " double" << std::endl;
    for (const auto& point : fluidPoints) {
        vtkFile << point.x << " " << point.y << " " << point.z << std::endl;
    }

    vtkFile << "CELLS " << fluidPoints.size() << " " << fluidPoints.size() * 2 << std::endl;
    for (size_t i = 0; i < fluidPoints.size(); ++i) {
        vtkFile << "1 " << i << std::endl;
    }

    vtkFile << "CELL_TYPES " << fluidPoints.size() << std::endl;
    for (size_t i = 0; i < fluidPoints.size(); ++i) {
        vtkFile << "1" << std::endl;
    }

    vtkFile << "POINT_DATA " << fluidPoints.size() << std::endl;

    vtkFile << "SCALARS i_index int 1" << std::endl;
    vtkFile << "LOOKUP_TABLE default" << std::endl;
    for (const auto& index : indicesI) {
        vtkFile << index << std::endl;
    }

    vtkFile << "SCALARS j_index int 1" << std::endl;
    vtkFile << "LOOKUP_TABLE default" << std::endl;
    for (const auto& index : indicesJ) {
        vtkFile << index << std::endl;
    }

    vtkFile << "SCALARS k_index int 1" << std::endl;
    vtkFile << "LOOKUP_TABLE default" << std::endl;
    for (const auto& index : indicesK) {
        vtkFile << index << std::endl;
    }

    vtkFile.close();
    std::cout << "VTK file " << filepath << " written successfully." << std::endl;
}


void Geometry::writeLabeledPointsToFile(const std::string& filename) {

#ifdef _WIN32
    std::filesystem::create_directory("labeledPoints");
#else 
    mkdir("labeledPoints", 0777);
#endif

    std::string filepath = "labeledPoints/labeled_points_" + std::to_string(myrank) + "_" + filename + ".vtk";

    std::ofstream vtkFile(filepath);

    if (!vtkFile.is_open()) {
        std::cerr << "Error: Unable to open VTK file for writing." << std::endl;
        return;
    }

    vtkFile << "# vtk DataFile Version 3.0" << std::endl;
    vtkFile << "Labeled Points Data" << std::endl;
    vtkFile << "ASCII" << std::endl;
    vtkFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

    std::vector<Point> fluidPoints;
    std::vector<int> indicesI, indicesJ, indicesK;

    for (size_t i = 0; i < pointsBox.size(); ++i) {
        for (size_t j = 0; j < pointsBox[i].size(); ++j) {
            for (size_t k = 0; k < pointsBox[i][j].size(); ++k) {
                const auto& point = pointsBox[i][j][k];
                if (point->PointBoxType == PointBoxType::real && point->NodeType == NodeType::Fluid) {
                    fluidPoints.push_back(*point);

                    std::vector<int> l2g = Local2Global(i, j, k);
                    indicesI.push_back(static_cast<int>(l2g[0]));
                    indicesJ.push_back(static_cast<int>(l2g[1]));
                    indicesK.push_back(static_cast<int>(l2g[2]));
                }
            }
        }
    }

    vtkFile << "POINTS " << fluidPoints.size() << " double" << std::endl;
    for (const auto& point : fluidPoints) {
        vtkFile << point.x << " " << point.y << " " << point.z << std::endl;
    }

    vtkFile << "CELLS " << fluidPoints.size() << " " << fluidPoints.size() * 2 << std::endl;
    for (size_t i = 0; i < fluidPoints.size(); ++i) {
        vtkFile << "1 " << i << std::endl;
    }

    vtkFile << "CELL_TYPES " << fluidPoints.size() << std::endl;
    for (size_t i = 0; i < fluidPoints.size(); ++i) {
        vtkFile << "1" << std::endl;
    }

    vtkFile << "POINT_DATA " << fluidPoints.size() << std::endl;

    vtkFile << "SCALARS i_index int 1" << std::endl;
    vtkFile << "LOOKUP_TABLE default" << std::endl;
    for (const auto& index : indicesI) {
        vtkFile << index << std::endl;
    }

    vtkFile << "SCALARS j_index int 1" << std::endl;
    vtkFile << "LOOKUP_TABLE default" << std::endl;
    for (const auto& index : indicesJ) {
        vtkFile << index << std::endl;
    }

    vtkFile << "SCALARS k_index int 1" << std::endl;
    vtkFile << "LOOKUP_TABLE default" << std::endl;
    for (const auto& index : indicesK) {
        vtkFile << index << std::endl;
    }

    vtkFile.close();
    std::cout << "VTK file " << filepath << " written successfully." << std::endl;
}


void Geometry::writePointsToFile(const std::vector<Point>& points, const std::string& foldername, const std::string& filename) {
#ifdef _WIN32
    std::filesystem::create_directory(foldername);
#else 
    mkdir(foldername.c_str(), 0777);
#endif

    std::string filepath = foldername + "/" + filename + "_" + std::to_string(myrank) + ".vtk";

    std::ofstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filepath << std::endl;
        return;
    }

    // Write the header information
    file << "# vtk DataFile Version 2.0" << std::endl;
    file << "Points Inside Complex Geometry" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;

    // Write the POINTS section
    int numPoints = points.size();
    file << "POINTS " << numPoints << " double" << std::endl;
    for (const auto& point : points) {
        file << point.x << " " << point.y << " " << point.z << std::endl;
    }

    // Write the CELLS section
    file << "CELLS " << numPoints << " " << numPoints * 2 << std::endl;
    for (int i = 0; i < numPoints; ++i) {
        file << "1 " << i << std::endl;
    }

    // Write the CELL_TYPES section
    file << "CELL_TYPES " << numPoints << std::endl;
    for (int i = 0; i < numPoints; ++i) {
        file << "1" << std::endl;
    }

    file.close();
}


