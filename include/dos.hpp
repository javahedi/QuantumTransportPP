#pragma once
#include "hamiltonian.hpp"
#include "mesh.hpp"
#include <vector>

/// @brief Density of States (DOS) calculator for solid-state systems
class DOS {
public:

    /// @brief Construct a DOS calculator with given Hamiltonian and mesh
    DOS(const Hamiltonian& H, const Mesh& mesh);
    // Compute total DOS over the specified energy window
    // energy_min, energy_max define the range
    // n_bins: number of energy bins
    // sigma: Gaussian smearing width (in same units as eigenvalues)

    /// @param Emin Minimum energy for DOS calculation
    /// @param Emax Maximum energy for DOS calculation
    /// @param nBins Number of energy bins for histogram
    /// @param sigma Gaussian smearing width (in same units as eigenvalues)
    /// @return Vector of DOS values for each energy bin
    /// @details The DOS is computed by integrating the density of states over the specified energy range.
    std::vector<double> computeDOS(double Emin, double Emax, int nBins, double sigma);
    /// @brief Compute projected DOS for a specific orbital index
    std::vector<double> computeProjectedDOS(double E_min, double E_max, int N_bins, double eta, int orbital_index);
    /// @brief Get the energy grid used for the DOS histogram
    const std::vector<double>& getEnergyGrid() const; 
   
private:
    const Hamiltonian& H;
    const Mesh& mesh;
    std::vector<double> energy_grid_; // Energies for the histogram

};