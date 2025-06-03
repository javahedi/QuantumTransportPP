#include "haldane.hpp"       // Your sample Hamiltonian (Haldane model)
#include "mesh.hpp"
#include "boltzmann.hpp"
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
    double tau = 1.0;
    double T = 0.01;           // Unitless (or Kelvin if enabled)
    double Ef = 0.0;           // Fermi level
    Eigen::Vector3d Efield(1.0, 0.0, 0.0); // E field in x direction

    // Construct Boltzmann solver
    BoltzmannSolver solver(H, mesh, tau,
                           false,  // temperature_in_kelvin
                           1.0);   // energy_scale

    Eigen::Matrix3d sigma = solver.conductivity(Ef, T, Efield);

    std::cout << "Conductivity tensor (σ_ij):\n" << sigma << std::endl;

    return 0;
}
