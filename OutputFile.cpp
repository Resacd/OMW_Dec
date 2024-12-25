#include "OutputFile.h"

#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkCellArray.h>
#include <vtkVertex.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkPolyDataWriter.h>


#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkDelaunay3D.h>
#include <vtkGeometryFilter.h>
#include <vtkPointData.h>
#include <vtkContourFilter.h>
#include <vtkSurfaceReconstructionFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPCANormalEstimation.h>

#include <vtkFillHolesFilter.h>
#include <vtkTriangle.h>



OutputFile::OutputFile(int outputIdx)
    : m_outputIdx(outputIdx)
{
}

OutputFile::~OutputFile()
{
}



void OutputFile::WriteVTK_LocalData(SubDomain& subdomain)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    std::string folder_name = "Fluid_Variables";

#ifdef _WIN32
    std::filesystem::create_directory(folder_name);
#else 
    mkdir(folder_name.c_str(), 0777);
#endif

    std::string file_name = folder_name + "/" + std::to_string(m_outputIdx) +
        "_rank_" + std::to_string(rank) + ".vtu";  

    vtkSmartPointer<vtkUnstructuredGrid> uGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();

    vtkSmartPointer<vtkDoubleArray> velocities = vtkSmartPointer<vtkDoubleArray>::New();
    velocities->SetNumberOfComponents(3);
    velocities->SetName("Velocity");
    vtkSmartPointer<vtkDoubleArray> pressures = vtkSmartPointer<vtkDoubleArray>::New();
    pressures->SetName("Pressure");
    vtkSmartPointer<vtkDoubleArray> velMags = vtkSmartPointer<vtkDoubleArray>::New();
    velMags->SetName("VelocityMagnitude");
    vtkSmartPointer<vtkDoubleArray> Density = vtkSmartPointer<vtkDoubleArray>::New();
    Density->SetName("Density");

    int pointId = 0;
    bool dataAdded = false;

    for (int z = 0; z < subdomain.my_Nz; ++z) {
        for (int y = 0; y < subdomain.my_Ny; ++y) {
            for (int x = 0; x < subdomain.my_Nx; ++x) {
                auto& node = subdomain.m_lattice[x][y][z];
                if (node != nullptr && node->NodeType_ != NodeType::Halo) {
                    double velMag = sqrt(node->velocity[0] * node->velocity[0] +
                        node->velocity[1] * node->velocity[1] +
                        node->velocity[2] * node->velocity[2]);

                   // std::vector<int> globalIndices = subdomain.Local2Global(node->x_position, node->y_position, node->z_position);
                   // points->InsertNextPoint(globalIndices[0], globalIndices[1], globalIndices[2]);
                    velocities->InsertNextTuple(node->velocity.data());
                    pressures->InsertNextValue(node->Pressure(subdomain.m_velSet));
                    velMags->InsertNextValue(velMag);
                    Density->InsertNextValue(node->density);

                    cells->InsertNextCell(1);
                    cells->InsertCellPoint(pointId);
                    ++pointId;

                    dataAdded = true;
                }
            }
        }
    }

    if (dataAdded) {
        uGrid->SetPoints(points);
        uGrid->SetCells(VTK_VERTEX, cells);

        // Add arrays to the point data
        uGrid->GetPointData()->AddArray(velocities);
        uGrid->GetPointData()->AddArray(pressures);
        uGrid->GetPointData()->AddArray(velMags);
        uGrid->GetPointData()->AddArray(Density);
    }

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(file_name.c_str());
    writer->SetInputData(uGrid);
    writer->Write();
}




void OutputFile::WriteVTK_Collective(SubDomain& subdomain)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    std::string folder_name = "Fluid_Variables";

#ifdef _WIN32
    std::filesystem::create_directory(folder_name);
#else
    mkdir(folder_name.c_str(), 0777);
#endif

    std::stringstream file_name_stream;
    file_name_stream << folder_name << "/" << std::to_string(m_outputIdx) << ".vtu";
    std::string file_name = file_name_stream.str();

    size_t local_size = 0;
    for (int z = 0; z < subdomain.my_Nz; ++z) {
        for (int y = 0; y < subdomain.my_Ny; ++y) {
            for (int x = 0; x < subdomain.my_Nx; ++x) {
                auto& node = subdomain.m_lattice[x][y][z];
                if (node != nullptr && node->NodeType_ != NodeType::Halo) {
                    local_size++;
                }
            }
        }
    }

    std::vector<double> local_data_send(local_size * 8);
    size_t index = 0;

    for (int z = 0; z < subdomain.my_Nz; ++z) {
        for (int y = 0; y < subdomain.my_Ny; ++y) {
            for (int x = 0; x < subdomain.my_Nx; ++x) {
                auto& node = subdomain.m_lattice[x][y][z];
                if (node != nullptr && node->NodeType_ != NodeType::Halo) {
                    double velMag = sqrt(node->velocity[0] * node->velocity[0] +
                        node->velocity[1] * node->velocity[1] +
                        node->velocity[2] * node->velocity[2]);
                    std::vector<int> globInD = subdomain.Geometry_.Local2Global(x, y, z);
                    local_data_send[index++] = static_cast<double>(globInD[0]);
                    local_data_send[index++] = static_cast<double>(globInD[1]);
                    local_data_send[index++] = static_cast<double>(globInD[2]);
                    local_data_send[index++] = node->velocity[0];
                    local_data_send[index++] = node->velocity[1];
                    local_data_send[index++] = node->velocity[2];
                    local_data_send[index++] = node->density;
                    local_data_send[index++] = velMag;
                }
            }
        }
    }

    std::vector<int> local_sizes(size);
    MPI_Gather(&local_size, 1, MPI_INT, local_sizes.data(), 1, MPI_INT, 0, comm);

    std::vector<int> displacements(size, 0);
    std::vector<int> recv_counts(size);
    if (rank == 0) {
        for (int i = 1; i < size; ++i) {
            displacements[i] = displacements[i - 1] + local_sizes[i - 1] * 8;
        }

        for (int i = 0; i < size; ++i) {
            recv_counts[i] = local_sizes[i] * 8;
        }
    }

    int total_size = std::accumulate(local_sizes.begin(), local_sizes.end(), 0) * 8;

    std::vector<double> all_data(total_size);

    MPI_Gatherv(local_data_send.data(), local_size * 8, MPI_DOUBLE,
        all_data.data(), recv_counts.data(), displacements.data(), MPI_DOUBLE, 0, comm);

    if (rank == 0) {
        vtkSmartPointer<vtkUnstructuredGrid> uGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
        vtkSmartPointer<vtkDoubleArray> velocityArray = vtkSmartPointer<vtkDoubleArray>::New();
        vtkSmartPointer<vtkDoubleArray> densityArray = vtkSmartPointer<vtkDoubleArray>::New();
        vtkSmartPointer<vtkDoubleArray> velMagArray = vtkSmartPointer<vtkDoubleArray>::New();

        velocityArray->SetNumberOfComponents(3);
        velocityArray->SetName("Velocity");

        densityArray->SetNumberOfComponents(1);
        densityArray->SetName("Density");

        velMagArray->SetNumberOfComponents(1);
        velMagArray->SetName("VelocityMagnitude");

        std::map<std::vector<int>, vtkIdType> pointIDMap;
        for (size_t i = 0; i < total_size / 8; ++i) {
            size_t dataIndex = i * 8;

            double x = all_data[dataIndex];
            double y = all_data[dataIndex + 1];
            double z = all_data[dataIndex + 2];

            vtkIdType pid = points->InsertNextPoint(x, y, z);
            std::vector<int> globalIndex = { static_cast<int>(x), static_cast<int>(y), static_cast<int>(z) };
            pointIDMap[globalIndex] = pid;

            double velocity[3] = { all_data[dataIndex + 3], all_data[dataIndex + 4], all_data[dataIndex + 5] };
            velocityArray->InsertNextTuple(velocity);

            double density = all_data[dataIndex + 6];
            densityArray->InsertNextValue(density);

            double velMag = all_data[dataIndex + 7];
            velMagArray->InsertNextValue(velMag);
        }

        int global_Nx = subdomain.Geometry_.nx;
        int global_Ny = subdomain.Geometry_.ny;
        int global_Nz = subdomain.Geometry_.nz;

        for (int z = 0; z < global_Nz - 1; ++z) {
            for (int y = 0; y < global_Ny - 1; ++y) {
                for (int x = 0; x < global_Nx - 1; ++x) {
                    std::vector<vtkIdType> cellPoints(8);

                    std::vector<int> indices[8] = {
                        {x, y, z},
                        {x + 1, y, z},
                        {x + 1, y + 1, z},
                        {x, y + 1, z},
                        {x, y, z + 1},
                        {x + 1, y, z + 1},
                        {x + 1, y + 1, z + 1},
                        {x, y + 1, z + 1}
                    };

                    bool validCell = true;
                    for (int i = 0; i < 8; ++i) {

                        if (pointIDMap.find(indices[i]) == pointIDMap.end()) {
                            validCell = false;
                            break;
                        }
                        cellPoints[i] = pointIDMap[indices[i]];
                    }

                    if (validCell) {
                        vtkSmartPointer<vtkHexahedron> hex = vtkSmartPointer<vtkHexahedron>::New();
                        for (int i = 0; i < 8; ++i) {
                            hex->GetPointIds()->SetId(i, cellPoints[i]);
                        }
                        uGrid->InsertNextCell(hex->GetCellType(), hex->GetPointIds());
                    }
                }
            }
        }



        uGrid->SetPoints(points);
        uGrid->GetPointData()->AddArray(velocityArray);
        uGrid->GetPointData()->AddArray(densityArray);
        uGrid->GetPointData()->AddArray(velMagArray);

        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(file_name.c_str());
        writer->SetInputData(uGrid);
        writer->Write();
    }
}


void OutputFile::WriteSurfaceAttribute(SubDomain& subdomain, int iter) {
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    std::string folder_name = "FlowHemodynamics";
    std::string file_name_prefix = "_FlowHemodynamics.csv";

#ifdef _WIN32
    std::filesystem::create_directory(folder_name);
#else 
    mkdir(folder_name.c_str(), 0777);
#endif

    size_t local_size = subdomain.nextToTheWall.size();
    std::vector<double> local_data_send(local_size * 6);  // 6 for X, Y, Z, WSS, OSI, TAWSS

    size_t index = 0;
    for (const auto& node : subdomain.nextToTheWall) {
        std::vector<int> ret = subdomain.Geometry_.RetrieveLocalInd(node->x_position, node->y_position, node->z_position);
        double x = subdomain.Geometry_.pointsBox[ret[0]][ret[1]][ret[2]]->x;
        double y = subdomain.Geometry_.pointsBox[ret[0]][ret[1]][ret[2]]->y;
        double z = subdomain.Geometry_.pointsBox[ret[0]][ret[1]][ret[2]]->z;

        local_data_send[index++] = x;
        local_data_send[index++] = y;
        local_data_send[index++] = z;
        local_data_send[index++] = node->wall_shear_stress;  // WSS
        local_data_send[index++] = node->osi;                // OSI
        local_data_send[index++] = node->tawss;              // TAWSS
    }

    std::vector<int> local_sizes(size);
    MPI_Gather(&local_size, 1, MPI_INT, local_sizes.data(), 1, MPI_INT, 0, comm);

    std::vector<int> displacements(size, 0);
    std::vector<int> recv_counts(size);
    if (rank == 0) {
        for (int i = 1; i < size; ++i) {
            displacements[i] = displacements[i - 1] + local_sizes[i - 1] * 6;
        }

        for (int i = 0; i < size; ++i) {
            recv_counts[i] = local_sizes[i] * 6;
        }
    }

    int total_size = std::accumulate(local_sizes.begin(), local_sizes.end(), 0) * 6;

    std::vector<double> all_data(total_size);

    MPI_Gatherv(local_data_send.data(), local_size * 6, MPI_DOUBLE,
        all_data.data(), recv_counts.data(), displacements.data(), MPI_DOUBLE, 0, comm);

    if (rank == 0) {
        std::string file_name = folder_name + "/" + std::to_string(iter) + file_name_prefix;
        std::ofstream outfile(file_name);

        if (outfile.is_open()) {
            outfile << "X,Y,Z,WSS,OSI,TAWSS\n";  
            for (size_t i = 0; i < total_size; i += 6) {
                outfile << all_data[i] << "," << all_data[i + 1] << "," << all_data[i + 2] << ","
                        << all_data[i + 3] << "," << all_data[i + 4] << "," << all_data[i + 5] << "\n";
            }
            outfile.close();
        }
        else {
            std::cerr << "Unable to open file " << file_name << std::endl;
        }
    }
}

