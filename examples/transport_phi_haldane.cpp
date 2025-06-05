#include "haldane.hpp"
#include "mesh.hpp"
#include "kubo.hpp"
#include "boltzmann.hpp"
#include "logInfo.hpp" 

int main() {
    // Parameters for the Haldane model
    const double t1 = 1.0;
    const double t2 = 0.1;
    const double M = 0.0;
    
    // 
    const double Ef = 1.0;
    const double eta = 1e-2;
    const double tau = 5.0;
    const double T = 0.01; // Temperature in Kelvin

    const bool temp_in_K = false;
    const double E_scale = t1; // Energy scale for the model

    const Eigen::Vector3d Efield(1.0, 0.0, 0.0);  // x-direction

    Mesh mesh(40, 40, 1, 1.0);

    std::ofstream out("transport_vs_phi.csv");
    out << "phi,kubo_xx,kubo_xy,boltzmann_xx,boltzmann_xy\n";

    logInfo("Starting phi loop from " + std::to_string(-M_PI) + " to " + std::to_string(M_PI));

    for (double phi = -M_PI; phi <= M_PI; phi += 0.1) {
        logInfo("Processing phi = " + std::to_string(phi));
        HaldaneModel H(t1, t2, M, phi);

        KuboSolver kubo(H, mesh, eta, temp_in_K, E_scale);
        Eigen::Matrix3d sig_kubo = kubo.conductivity(Ef, T);

        BoltzmannSolver bolt(H, mesh, tau, temp_in_K, E_scale);
        Eigen::Matrix3d sig_bolt = bolt.conductivity(Ef, T, Efield);

        out << phi << ","
            << sig_kubo(0, 0) << "," << sig_kubo(0, 1) << ","
            << sig_bolt(0, 0) << "," << sig_bolt(0, 1) << "\n";
    }

    out.close();
    logSuccess("Transport calculations completed successfully.");
    return 0;
}
