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


void Writer::writeLabeledPointsToFile(const std::vector<std::vector<std::vector<Point>>>& pointsBox, const std::string& filename) {
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
                if (point.Type == NodeType::Fluid) {
                    fluidPoints.push_back(point);
                }
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

void Writer::writeLineToVTK(const Line& line, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    // Write the header information
    file << "# vtk DataFile Version 2.0" << std::endl;
    file << "Line for Visualization" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;

    // Write the POINTS section
    file << "POINTS 2 double" << std::endl;
    file << line.start.x << " " << line.start.y << " " << line.start.z << std::endl;
    file << line.end.x << " " << line.end.y << " " << line.end.z << std::endl;

    // Write the CELLS section for the line
    file << "CELLS 1 3" << std::endl;
    file << "2 0 1" << std::endl;

    // Write the CELL_TYPES section (use VTK_LINE)
    file << "CELL_TYPES 1" << std::endl;
    file << "3" << std::endl;

    file.close();
}
