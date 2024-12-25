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
#include "Utility.h"
#include "GeometryReader.h"
#include "Input.h"

int main(int argc, char* argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    int myrank, num_procs;
    Utility Utility;
    Reader Reader;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    if (myrank == 0) {
        Print_disclaimer();

        omp_set_num_threads(omp_get_max_threads());
        int numThreads = omp_get_max_threads();
        omp_set_num_threads(numThreads);
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
        chronoThread.detach();
    }


    double dx = 2e-5;
    double diameter = 1e-2;
    double dt = 1.8e-5;
    double rho = 1025.0;
    double nu_blood = 3.8e-6;

    double u_LB = u_in * (dt / dx);
    double nu_LB = nu_blood * (dt / (std::pow(dx, 2)));
   // double tau = 3 * nu_LB + 0.5;
    double tau=0.68;
    //double tau = 10;

    double Reynolds = u_in * diameter / nu_blood;

if(myrank==0){
   // std::cout << "Relaxation Time is:" << " " << tau << std::endl;
   // std::cout << "Reylonds is: " << " " << Reynolds << std::endl;
}
    std::shared_ptr<LatticeModel> velSet = std::make_shared<LatticeModel>();
    Geometry Geometry("OMW2011.stl",velSet, spacing, 'z');
    SubDomain SubDomain(Geometry);
    Simulation(SubDomain, tau, pressure,dt);

    MPI_Finalize();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    if (myrank == 0) {
        std::cout << "\rElapsed time: " << std::fixed << std::setprecision(2)
            << std::setw(8) << elapsed.count() << " seconds." << std::endl;
    }

    return 0;
}
