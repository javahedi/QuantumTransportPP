#pragma once
#include "hamiltonian.hpp"
#include "mesh.hpp"
#include <vector>

class DOS {
public:
    DOS(const Hamiltonian& H, const Mesh& mesh);
    // Compute total DOS over the specified energy window
    // energy_min, energy_max define the range
    // n_bins: number of energy bins
    // sigma: Gaussian smearing width (in same units as eigenvalues)

    std::vector<double> computeDOS(double Emin, double Emax, int nBins, double sigma);
    std::vector<double> computeProjectedDOS(double E_min, double E_max, int N_bins, double eta, int orbital_index);

    const std::vector<double>& getEnergyGrid() const; 
   
private:
    const Hamiltonian& H;
    const Mesh& mesh;
    std::vector<double> energy_grid_; // Energies for the histogram

};