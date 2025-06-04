#include "haldane.hpp"       // Your sample Hamiltonian (Haldane model)
#include "mesh.hpp"
#include "dos.hpp"
#include <iostream>
#include <Eigen/Dense>

/*

DOS(E) = \sum_{n,k}\delta(E-\epsilon_{nk})
PDOS_i(E) = \sum_{n,k}\delta(E-\epsilon_{nk})|<i|\psi_{nk}>|^2

*/
int main() {
    // Set up k-mesh
    size_t N = 30; // Grid resolution
    Mesh mesh(N, N); // 2D mesh

    // Haldane model parameters
    HaldaneModel H;
    H.t1 = 1.0;
    H.t2 = 0.1;
    H.phi = M_PI / 2.0;
    H.M = 0.2; 

    // Relaxation time and temperature
    double Ef = 0.0;            // Fermi energy (unitless)
    double eta = 1e-2;          // Broadening

    // Construct Boltzmann solver
    DOS dos(H, mesh);   
    std::vector<double> dos_vals = dos.computeDOS(Ef - 4.0, Ef + 4.0, 400, eta);
    std::vector<double> energy_grid = dos.getEnergyGrid();

    // Output the results
    std::cout << "# Energy\tDOS\n";
    for (size_t i = 0; i < dos_vals.size(); ++i) {
        std::cout << energy_grid[i] << "\t" << dos_vals[i] << "\n";
    }

    return 0;
}
