/** \file main.cpp
 * \brief Computes conductivity for an altermagnet model using Kubo and Boltzmann solvers.
 * \details Iterates over J values (0.0 to 0.10, step 0.01) and lambda values, outputting
 *          xx and xy conductivity components to separate CSV files for each lambda
 *          (e.g., data/altermagnet_conductivity_lambda_0.0.csv).
 * \author [Your Name]
 * \date June 2025
 */

#include "mesh.hpp"
#include "altermagnet.hpp"
#include "geometry.hpp"
#include "logInfo.hpp"
#include "boltzmann.hpp"
#include "kubo.hpp"
#include <vector>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <omp.h>

using namespace QuantumTransport;


int main(int argc, char* argv[])  {

    //omp_set_num_threads(4);  // Set to desired number of threads
    if (argc > 1) {
        int threads = std::atoi(argv[1]);
        omp_set_num_threads(threads);
    }

    // Initialize parameters
    std::vector<double> Jlist;

    for (double j = 0.0; j <= 1.0 ; j += 0.025) { // 1e-6 ensures 0.10 is included
        Jlist.push_back(j);
    }

    std::vector<double> lambdaList = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5}; ///< Lambda values for iteration
    size_t N = 250; ///< Grid resolution
    Mesh mesh(N, N, 1, M_PI); ///< 2D mesh for calculations

    logInfo("Starting altermagnet conductivity calculation...");

    // Kubo and Boltzmann solver parameters
    double Ef  = 0.5; ///< Fermi energy (unitless)
    double T   = 0.02; ///< Temperature (unitless)
    double eta = 1e-2; ///< Broadening for Kubo solver
    double tau = 100.0; ///< Relaxation time for Boltzmann solver
    Eigen::Vector3d Efield(1.0, 0.0, 0.0); ///< Electric field in x direction
    Eigen::Vector3d gradT(1.0, 0.0, 0.0); ///< Temperature gradient in x direction
    Eigen::Vector3d Bfield(0, 0, 1.0);   ///< Magnetic field in z direction

    #pragma omp parallel for
    for (const auto& lambda : lambdaList) {

        // Create unique filename for each lambda
        std::stringstream filename;
        filename << std::fixed << std::setprecision(1) << "altermagnet_conductivity_EF0.5_lambda_" << lambda << ".csv";
        logInfo("Writing results to " + filename.str() + "...");

        // Open output file for this lambda
       std::ofstream out(filename.str());
        if (!out) {  // More idiomatic check
            logError("Failed to open " + filename.str() + ": " + strerror(errno));
            continue;
        }
        // Write header (flushed once)
        {
            out << std::fixed << std::setprecision(6);
            out << "J,sigma_xx,sigma_xy,alpha_xx,alpha_xy\n";
            out.flush();
        }

        
        for (const auto& j : Jlist) {

            AltermagnetModel altermagnet(1.0, j, lambda); ///< Altermagnet model with updated J and lambda

            // Kubo solver
            //KuboSolver kubo(altermagnet, mesh, eta, false, 1.0);
            //auto [sigmaK, alphaK, kappa] = kubo.computeTransportTensors(Ef, T);

            // Boltzmann solver
            BoltzmannSolver boltzmann(altermagnet, mesh, tau, false, 1.0);
            auto [sigma, alpha] = boltzmann.computeTransportTensors(Ef, T, gradT, Efield, Bfield);

            // compute Seebeck sigma.inverse() * alpha
            // Eigen::Matrix3d Seebeck = sigma.inverse() * alpha; 
            
            // Write results
            out << j << ","
                << sigma(0, 0) << "," << sigma(0, 1) << ","
                << alpha(0, 0) << "," << alpha(0, 1) << "\n";
        }

        out.close();
    }

    logSuccess("Altermagnet conductivity calculation completed successfully.");
    return 0;
}