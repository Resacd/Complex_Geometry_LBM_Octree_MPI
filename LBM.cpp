#include <iostream>
#include "LatticeModel.h"
#include "Simulation.h"
#include "OutputFile.h"
#include <mpi.h>
#include "Subdomain.h"
#include "Subdomain.h"
#include "Input.h"
#include "ComplexGeometry.h"
#include "disclaimer.h"
#include <iostream>
#include <chrono>
#include <thread>
#include <iomanip>
#include <omp.h>


int main(int argc, char* argv[])
{

    auto start = std::chrono::high_resolution_clock::now();
    int myrank, num_procs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    if (myrank==0){ Print_disclaimer(); }
    omp_set_num_threads(omp_get_max_threads());
    int numThreads = omp_get_max_threads();
    omp_set_num_threads(numThreads);
    std::cout << "Number of OpenMP threads: " << numThreads << std::endl;

    bool stopChrono = false;
    std::thread chronoThread([&]() {
        while (!stopChrono) {
            auto now = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = now - start;

            std::cout << "\rElapsed time: " << std::fixed << std::setprecision(2)
                << std::setw(8) << elapsed.count() << " seconds" << std::flush;

            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        });

    std::shared_ptr<LatticeModel> velSet = std::make_shared<LatticeModel>();



    Geometry Geometry("CylinderASCII.stl", InOutCoord, spacing, mdepth);
    //Geometry Geometry("bifurcation.stl",InOutCoord, spacing,mdepth);
    //SubDomain SubDomain(velSet, Geometry);
    //Simulation(SubDomain, tau, inVel, pressure);
    stopChrono = true;
    chronoThread.join();

    MPI_Finalize();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "\rElapsed time: " << std::fixed << std::setprecision(2)
        << std::setw(8) << elapsed.count() << " seconds." << std::endl;

    return 0;
}



