#include "OutputFile.h"
#include "Struct.h"
#include <HaloLayer.h>

OutputFile::OutputFile(int outputIdx)
    : m_outputIdx(outputIdx)
{
}

OutputFile::~OutputFile()
{
}


void OutputFile::WriteVTK(std::vector<std::vector<std::vector<std::shared_ptr<Node>>>>& lattice,
    std::shared_ptr<VelocitySet> velSet, int size_x, int size_y, int size_z,
    int iter)
{

    std::string vtkFileName = "output_" + std::to_string(iter) + ".vtk";

    std::ofstream vtkFile(vtkFileName);

    if (!vtkFile.is_open()) {
        std::cerr << "Error: Could not open VTK file for writing." << std::endl;
        return;
    }

    vtkFile << "# vtk DataFile Version 3.0" << std::endl;
    vtkFile << "Lattice Boltzmann Simulation Data" << std::endl;
    vtkFile << "ASCII" << std::endl;
    vtkFile << "DATASET STRUCTURED_POINTS" << std::endl;
    vtkFile << "DIMENSIONS " << size_x << " " << size_y << " " << size_z << std::endl;
    vtkFile << "ORIGIN 0 0 0" << std::endl;
    vtkFile << "SPACING 1 1 1" << std::endl;
    vtkFile << "POINT_DATA " << size_x * size_y * size_z << std::endl;

    vtkFile << "SCALARS VelocityX double" << std::endl;
    vtkFile << "LOOKUP_TABLE default" << std::endl;

    for (int z = 0; z < size_z; ++z) {
        for (int y = 0; y < size_y; ++y) {
            for (int x = 0; x < size_x; ++x) {
                auto& node = lattice[x][y][z];
                if (node != nullptr) {
                    vtkFile << node->Velocity(velSet)[0] << std::endl;
                }
                else {
                    vtkFile << "0.0" << std::endl;
                }
            }
        }
    }

    vtkFile << "SCALARS VelocityY double" << std::endl;
    vtkFile << "LOOKUP_TABLE default" << std::endl;

    for (int z = 0; z < size_z; ++z) {
        for (int y = 0; y < size_y; ++y) {
            for (int x = 0; x < size_x; ++x) {
                auto& node = lattice[x][y][z];
                if (node != nullptr) {
                    vtkFile << node->Velocity(velSet)[1] << std::endl;
                }
                else {
                    vtkFile << "0.0" << std::endl;
                }
            }
        }
    }

    vtkFile << "SCALARS VelocityZ double" << std::endl;
    vtkFile << "LOOKUP_TABLE default" << std::endl;

    for (int z = 0; z < size_z; ++z) {
        for (int y = 0; y < size_y; ++y) {
            for (int x = 0; x < size_x; ++x) {
                auto& node = lattice[x][y][z];
                if (node != nullptr) {
                    vtkFile << node->Velocity(velSet)[2] << std::endl;
                }
                else {
                    vtkFile << "0.0" << std::endl;
                }
            }
        }
    }

    vtkFile << "SCALARS Pressure double" << std::endl;
    vtkFile << "LOOKUP_TABLE default" << std::endl;

    for (int z = 0; z < size_z; ++z) {
        for (int y = 0; y < size_y; ++y) {
            for (int x = 0; x < size_x; ++x) {
                auto& node = lattice[x][y][z];
                if (node != nullptr) {
                    vtkFile << node->Pressure(velSet) << std::endl;
                }
                else {
                    vtkFile << "0.0" << std::endl;
                }
            }
        }
    }

    vtkFile << "SCALARS VelocityMagnitude double" << std::endl;
    vtkFile << "LOOKUP_TABLE default" << std::endl;

    for (int z = 0; z < size_z; ++z) {
        for (int y = 0; y < size_y; ++y) {
            for (int x = 0; x < size_x; ++x) {
                auto& node = lattice[x][y][z];
                if (node != nullptr) {
                    double velMag = sqrt(node->Velocity(velSet)[0] * node->Velocity(velSet)[0] +
                        node->Velocity(velSet)[1] * node->Velocity(velSet)[1] +
                        node->Velocity(velSet)[2] * node->Velocity(velSet)[2]);
                    vtkFile << velMag << std::endl;
                }
                else {
                    vtkFile << "0.0" << std::endl; 
                }
            }
        }
    }

    vtkFile.close();

    std::cout << "VTK file '" << vtkFileName << "' written successfully." << std::endl;
}


void OutputFile::WriteCSV(std::vector<std::vector<std::vector<std::shared_ptr<Node>>>>& lattice,
    std::shared_ptr<VelocitySet> velSet, int start_x, int start_y, int start_z, int end_x, int end_y, int end_z)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    std::string file_name = std::to_string(m_outputIdx) + ".csv";
    std::ofstream velocityCSV;

    if (rank == 0) {
        velocityCSV.open(file_name, std::ios_base::app); 

        if (!velocityCSV.is_open()) {
            std::cerr << "Rank " << rank << ": Error opening file " << file_name << std::endl;
            MPI_Abort(comm, 1);
        }

        velocityCSV << std::setw(10) << "coordX" << ","
            << std::setw(10) << "coordY" << ","
            << std::setw(10) << "coordZ" << ","
            << std::setw(10) << "velX" << ","
            << std::setw(10) << "velY" << ","
            << std::setw(10) << "velZ" << ","
            << std::setw(10) << "pressure" << ","
            << std::setw(10) << "velMag" << std::endl;

        velocityCSV.close(); 
    }

    MPI_Barrier(comm);

    velocityCSV.open(file_name, std::ios_base::app); 

    if (!velocityCSV.is_open()) {
        std::cerr << "Rank " << rank << ": Error opening file " << file_name << std::endl;
        MPI_Abort(comm, 1);
    }

    std::vector<std::vector<double>> all_data;

    for (int z = start_z; z <= end_z; ++z) {
        for (int y = start_y; y <= end_y; ++y) {
            for (int x = start_x; x <= end_x; ++x) {
                int global_x = x;
                int global_y = y;
                int global_z = z;
                auto& node = lattice[x - start_x][y - start_y][z - start_z];
                if (node != nullptr && node->NodeType_!=NodeType::Halo) {
                    double velMag = sqrt(node->Velocity(velSet)[0] * node->Velocity(velSet)[0] +
                        node->Velocity(velSet)[1] * node->Velocity(velSet)[1] + node->Velocity(velSet)[2] *
                        node->Velocity(velSet)[2]);

                    // Store data in the vector
                    all_data.push_back({ static_cast<double>(global_x),
                        static_cast<double>(global_y),
                        static_cast<double>(global_z),
                        node->Velocity(velSet)[0],
                        node->Velocity(velSet)[1],
                        node->Velocity(velSet)[2],
                        node->Pressure(velSet),
                        velMag });
                }
            }
        }
    }

    std::sort(all_data.begin(), all_data.end(),
        [](const std::vector<double>& a, const std::vector<double>& b) {
            return std::tie(a[0], a[1], a[2]) < std::tie(b[0], b[1], b[2]);
        });

    MPI_Barrier(comm);

    for (const auto& data : all_data) {
        velocityCSV << std::fixed << std::setprecision(6);

        velocityCSV << std::setw(10) << static_cast<int>(data[0]) << ","
            << std::setw(10) << static_cast<int>(data[1]) << ","
            << std::setw(10) << static_cast<int>(data[2]) << ","
            << std::setw(10) << data[3] << ","
            << std::setw(10) << data[4] << ","
            << std::setw(10) << data[5] << ","
            << std::setw(10) << data[6] << ","
            << std::setw(10) << data[7] << std::endl;
    }

    MPI_Barrier(comm);

    if (rank == 0) {
        velocityCSV.close();
    }
}




