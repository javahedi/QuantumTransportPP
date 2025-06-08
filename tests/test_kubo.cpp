#include "haldane.hpp"       // Your sample Hamiltonian (Haldane model)
#include "mesh.hpp"
#include "kubo.hpp"
#include <iostream>
#include <Eigen/Dense>

int main() {
    // Set up k-mesh
    size_t N = 30; // Grid resolution
    Mesh mesh(N, N); // 2D mesh

    // Haldane model parameters
    double t1 = 1.0;
    double t2 = 0.1;
    double phi = M_PI / 2.0;

    HaldaneModel H(t1, t2, phi);

    // Relaxation time and temperature
    double Ef = 0.0;            // Fermi energy (unitless)
    double T = 0.01;            // Temperature (unitless)
    double eta = 1e-2;          // Broadening

    // Construct Boltzmann solver
    KuboSolver solver(H, mesh, eta,
                           false,  // temperature_in_kelvin
                           1.0);   // energy_scale

    auto [sigma, alpha, kappa] = solver.computeTransportTensors(Ef, T);

    std::cout << "Conductivity tensor (σ_ij):\n" << sigma << std::endl;

    return 0;
}
