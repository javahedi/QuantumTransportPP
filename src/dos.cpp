#include "dos.hpp"
#include <cmath>
#include <Eigen/Eigenvalues>
#include <iostream>

DOS::DOS(const Hamiltonian& H, const Mesh& mesh)
    : H(H), mesh(mesh) {}

std::vector<double> DOS::computeDOS(double energy_min, double energy_max, int n_bins, double sigma) {
    energy_grid_.resize(n_bins);
    std::vector<double> dos(n_bins, 0.0);

    double dE = (energy_max - energy_min) / n_bins;
    for (int i = 0; i < n_bins; ++i) {
        energy_grid_[i] = energy_min + (i + 0.5) * dE;
    }

    // Loop over all k-points
    for (const auto& k : mesh.getKPoints()) {
        Eigen::VectorXd evals;
        Eigen::MatrixXcd evecs;
        H.eigensystem(k, evals, evecs);

        // For each eigenvalue, add Gaussian smeared delta peak to DOS
        for (int n = 0; n < evals.size(); ++n) {
            double energy = evals[n];
            // Add contribution to all bins with significant weight (within ~3*sigma)
            int start_bin = std::max(0, int((energy - 3 * sigma - energy_min) / dE));
            int end_bin = std::min(n_bins - 1, int((energy + 3 * sigma - energy_min) / dE));

            for (int bin = start_bin; bin <= end_bin; ++bin) {
                double x = (energy_grid_[bin] - energy) / sigma;
                double weight = std::exp(-0.5 * x * x) / (sigma * std::sqrt(2 * M_PI));
                dos[bin] += weight;
            }
        }
    }

    // Normalize by number of k-points and bin width to get DOS per energy unit
    double norm = 1.0 / (mesh.size() * dE);
    for (auto& val : dos) {
        val *= norm;
    }

    return dos;
}


std::vector<double> DOS::computeProjectedDOS(double E_min, double E_max, int N_bins, double eta, int orbital_index) {
    std::vector<double> pdos(N_bins, 0.0);
    energy_grid_.resize(N_bins);
    double dE = (E_max - E_min) / N_bins;

    for (int i = 0; i < N_bins; ++i) {
        energy_grid_[i] = E_min + (i + 0.5) * dE;
    }

    Eigen::VectorXd evals;
    Eigen::MatrixXcd evecs;

    for (const auto& k : mesh.getKPoints()) {
        H.eigensystem(k, evals, evecs);
        for (int n = 0; n < evals.size(); ++n) {
            std::complex<double> amp = evecs.col(n)(orbital_index);
            double weight = std::norm(amp);  // |⟨i|ψ⟩|^2

            for (int i = 0; i < N_bins; ++i) {
                double x = (energy_grid_[i] - evals[n]) / eta;
                pdos[i] += weight / (M_PI * eta * (1 + x * x));
            }
        }
    }

    double norm = 1.0 / static_cast<double>(mesh.size());
    for (auto& d : pdos) {
        d *= norm;
    }

    return pdos;
}

const std::vector<double>& DOS::getEnergyGrid() const {
    return energy_grid_;
}
