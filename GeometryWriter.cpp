#include "GeometryWriter.h"


Writer::Writer() {


}

Writer::~Writer() {


}

void Writer::writePointoffToFile(const std::vector<Point>& points, const std::string& filename) {
    // Open the output file for writing
    std::ofstream outFile(filename);

    // Write the header information for the VTK file format
    outFile << "# vtk DataFile Version 3.0" << std::endl;
    outFile << "Points data" << std::endl;
    outFile << "ASCII" << std::endl;
    outFile << "DATASET POLYDATA" << std::endl;

    // Write the points data
    outFile << "POINTS " << points.size() << " double" << std::endl;
    for (const auto& point : points) {
        outFile << point.x << " " << point.y << " " << point.z << std::endl;
    }

    // Close the output file
    outFile.close();
}


void Writer::writePointsToFile(const std::vector<Point>& points, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
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



void Writer::writeBoxToFile(const std::vector<std::vector<std::vector<Point>>>& points, const std::string& filename) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Unable to open file " << filename << " for writing" << std::endl;
        return;
    }

    int nx = points.size();
    int ny = points[0].size();
    int nz = points[0][0].size();

    int numPoints = nx * ny * nz;

    file << "# vtk DataFile Version 3.0\n";
    file << "Points\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";
    file << "POINTS " << numPoints << " float\n";

    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                file << (points[i][j][k].x) << " " << (points[i][j][k].y) << " " << (points[i][j][k].z) << "\n";
            }
        }
    }


    file.close();
    std::cout << " " << "Points written to file " << filename << std::endl;
}

void Writer::writePointsToFile(const std::vector<std::vector<std::vector<Point>>>& pointsBox, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    // Write the header information
    file << "# vtk DataFile Version 2.0" << std::endl;
    file << "Points Inside Complex Geometry" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;

    // Collect the fluid points
    std::vector<Point> fluidPoints;
    for (const auto& plane : pointsBox) {
        for (const auto& row : plane) {
            for (const auto& point : row) {

                    fluidPoints.push_back(point);
            }
        }
    }

    // Write the POINTS section
    file << "POINTS " << fluidPoints.size() << " double" << std::endl;
    for (const auto& point : fluidPoints) {
        file << point.x << " " << point.y << " " << point.z << std::endl;
    }

    // Write the CELLS section
    int numCells = fluidPoints.size();
    file << "CELLS " << numCells << " " << numCells * 2 << std::endl;
    for (int i = 0; i < numCells; ++i) {
        file << "1 " << i << std::endl;
    }

    // Write the CELL_TYPES section
    file << "CELL_TYPES " << numCells << std::endl;
    for (int i = 0; i < numCells; ++i) {
        file << "1" << std::endl;
    }

    file.close();
}

void Writer::writeLabeledPointsToFile(std::vector< std::vector< std::vector<std::shared_ptr<Point>>>>& pointsBox, const std::string& filename) {
    std::ofstream vtkFile(filename);

    if (!vtkFile.is_open()) {
        std::cerr << "Error: Unable to open VTK file for writing." << std::endl;
        return;
    }

    vtkFile << "# vtk DataFile Version 3.0" << std::endl;
    vtkFile << "Fluid Points Data" << std::endl;
    vtkFile << "ASCII" << std::endl;
    vtkFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

    // Collect the fluid points and their indices
    std::vector<Point> fluidPoints;
    std::vector<int> indicesI, indicesJ, indicesK;

    for (size_t i = 0; i < pointsBox.size(); ++i) {
        for (size_t j = 0; j < pointsBox[i].size(); ++j) {
            for (size_t k = 0; k < pointsBox[i][j].size(); ++k) {
                const auto& point = pointsBox[i][j][k];
                if (point) {
                    fluidPoints.push_back(*point);
                    indicesI.push_back(i);
                    indicesJ.push_back(j);
                    indicesK.push_back(k);
                }
            }
        }
    }

    // Write the POINTS section
    vtkFile << "POINTS " << fluidPoints.size() << " double" << std::endl;
    for (const auto& point : fluidPoints) {
        vtkFile << point.x << " " << point.y << " " << point.z << std::endl;
    }

    // Write the CELLS section
    vtkFile << "CELLS " << fluidPoints.size() << " " << fluidPoints.size() * 2 << std::endl;
    for (size_t i = 0; i < fluidPoints.size(); ++i) {
        vtkFile << "1 " << i << std::endl;
    }

    // Write the CELL_TYPES section
    vtkFile << "CELL_TYPES " << fluidPoints.size() << std::endl;
    for (size_t i = 0; i < fluidPoints.size(); ++i) {
        vtkFile << "1" << std::endl;
    }

    // Add point data fields for i, j, k indices
    vtkFile << "POINT_DATA " << fluidPoints.size() << std::endl;

    // Writing i indices
    vtkFile << "SCALARS i_index int 1" << std::endl;
    vtkFile << "LOOKUP_TABLE default" << std::endl;
    for (const auto& index : indicesI) {
        vtkFile << index << std::endl;
    }

    // Writing j indices
    vtkFile << "SCALARS j_index int 1" << std::endl;
    vtkFile << "LOOKUP_TABLE default" << std::endl;
    for (const auto& index : indicesJ) {
        vtkFile << index << std::endl;
    }

    // Writing k indices
    vtkFile << "SCALARS k_index int 1" << std::endl;
    vtkFile << "LOOKUP_TABLE default" << std::endl;
    for (const auto& index : indicesK) {
        vtkFile << index << std::endl;
    }

    vtkFile.close();
    std::cout << "VTK file " << filename << " written successfully." << std::endl;
}


