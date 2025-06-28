/** \file altermagnet_conductivity.cpp
 * \brief Computes conductivity for an altermagnet model using Kubo and Boltzmann solvers.
 * \details Iterates over J values (0.0 to 0.10, step 0.01) and lambda values, outputting
 *          xx and xy conductivity components to separate CSV files for each lambda
 *          (e.g., data/altermagnet_conductivity_lambda_0.0.csv).
 * \author Javad Vahedi
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
#include <filesystem>

namespace fs = std::filesystem;
using namespace QuantumTransport;

int main(int argc, char* argv[]) {
    // Set OpenMP threads
    if (argc > 1) {
        int threads = std::atoi(argv[1]);
        omp_set_num_threads(threads);
    }

    // Initialize parameters
    std::vector<double> Jlist;
    for (double j = 0.0; j <= 1.0; j += 0.1) {
        Jlist.push_back(j);
    }

    std::vector<double> lambdaList = {0.0, 0.1, 0.2, 0.3};
    size_t N = 100;
    Mesh mesh(N, N, 1, M_PI);

    // Transport parameters
    double Ef = 0.5;
    double T = 0.02;
    double eta = 1e-2;
    double tau = 100.0;
    Eigen::Vector3d Efield(1.0, 0.0, 0.0);
    Eigen::Vector3d gradT(1.0, 0.0, 0.0);
    Eigen::Vector3d Bfield(0, 0, 1.0);

    // Create data directory (relative to executable location)
    std::string output_dir = "../data/";
    fs::create_directories(output_dir);

    logInfo("Starting altermagnet conductivity calculation...");

    #pragma omp parallel for
    for (const auto& lambda : lambdaList) {
        std::stringstream filename;
        filename << std::fixed << std::setprecision(1) 
                 << "altermagnet_conductivity_EF" << Ef << "_lambda_" << lambda << ".csv";
        
        std::string output_path = output_dir + filename.str();
        logInfo("Writing results to " + output_path);

        std::ofstream out(output_path);
        if (!out) {
            logError("Failed to open " + output_path + ": " + strerror(errno));
            continue;
        }

        // Write header
        out << std::fixed << std::setprecision(6);
        out << "J,sigma_xx,sigma_xy,alpha_xx,alpha_xy\n";
        out.flush();

        for (const auto& j : Jlist) {
            AltermagnetModel altermagnet(1.0, j, lambda);
            BoltzmannSolver boltzmann(altermagnet, mesh, tau, false, 1.0);
            auto [sigma, alpha] = boltzmann.computeTransportTensors(Ef, T, gradT, Efield, Bfield);

            out << j << ","
                << sigma(0, 0) << "," << sigma(0, 1) << ","
                << alpha(0, 0) << "," << alpha(0, 1) << "\n";
        }
        out.close();
    }

    logSuccess("Altermagnet conductivity calculation completed successfully.");
    return 0;
}