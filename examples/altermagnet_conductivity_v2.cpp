/** \file altermagnet_conductivity.cpp
 * \brief Computes conductivity for an altermagnet model using Kubo and Boltzmann solvers.
 * \details Reads parameters from params.ini file and outputs conductivity components to CSV files in the data/ directory.
 * \author Javad Vahedi
 * \date June 2025
 */

#include "mesh.hpp"
#include "altermagnet.hpp"
#include "geometry.hpp"
#include "logInfo.hpp"
#include "boltzmann.hpp"
#include "kubo.hpp"
#include "config_parser.hpp"
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

    try {
        // Load configuration
        ConfigParser config("../params.ini");
        
        // Get system parameters
        auto J_values = config.getRange("System", "J", {0.0, 0.5, 1.0});
        auto lambda_values = config.getRange("System", "lambda", {0.0, 0.1, 0.5});
        size_t mesh_points = config.get<size_t>("System", "mesh_points", 250);
        
        // Get transport parameters
        double Ef = config.get<double>("Transport", "Ef", 0.5);
        double temperature = config.get<double>("Transport", "temperature", 0.02);
        double eta = config.get<double>("Transport", "eta", 1e-2);
        double tau = config.get<double>("Transport", "tau", 100.0);

        // Field configurations
        Eigen::Vector3d Efield(1.0, 0.0, 0.0);  // Electric field in x direction
        Eigen::Vector3d gradT(1.0, 0.0, 0.0);    // Temperature gradient in x direction  
        Eigen::Vector3d Bfield(0, 0, 1.0);       // Magnetic field in z direction

        // Initialize mesh
        Mesh mesh(mesh_points, mesh_points, 1, M_PI);
        
        // Create data directory (relative to executable location)
        std::string output_dir = "../data/";
        fs::create_directories(output_dir);

        logInfo("Starting altermagnet conductivity calculation...");

        #pragma omp parallel for
        for (const auto& lambda : lambda_values) {
            // Create output filename
            std::stringstream filename;
            filename << std::fixed << std::setprecision(1) 
                    << "altermagnet_conductivity_EF" << Ef << "_lambda_" << lambda << ".csv";
            
            std::string output_path = output_dir + filename.str();
            logInfo("Writing results to " + output_path);

            // Open output file
            std::ofstream out(output_path);
            if (!out) {
                logError("Failed to open " + output_path + ": " + strerror(errno));
                continue;
            }

            // Write header
            out << std::fixed << std::setprecision(6);
            out << "J,sigma_xx,sigma_xy,alpha_xx,alpha_xy\n";
            out.flush();

            // Main computation loop
            for (const auto& j : J_values) {
                AltermagnetModel altermagnet(1.0, j, lambda);

                // Boltzmann solver
                BoltzmannSolver boltzmann(altermagnet, mesh, tau, false, 1.0);
                auto [sigma, alpha] = boltzmann.computeTransportTensors(Ef, temperature, gradT, Efield, Bfield);

                // Write results
                out << j << ","
                    << sigma(0, 0) << "," << sigma(0, 1) << ","
                    << alpha(0, 0) << "," << alpha(0, 1) << "\n";
            }
            out.close();
        }

        logSuccess("Calculation completed successfully. Output files saved in data/ directory.");
        return 0;

    } catch (const std::exception& e) {
        logError("Error: " + std::string(e.what()));
        return 1;
    }
}