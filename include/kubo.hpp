
#pragma once
#include "hamiltonian.hpp"
#include "mesh.hpp"
#include <Eigen/Dense>

/// @brief Kubo-Greenwood solver for calculating the conductivity tensor
class KuboSolver {
public:
    /// @brief Construct a KuboSolver with given Hamiltonian and mesh
    KuboSolver(const Hamiltonian& H, const Mesh& mesh, double eta = 1e-3,
               bool temperature_in_kelvin = false, double energy_scale = 1.0);

    /// @brief Compute conductivity, thermoelectric, and thermal conductivity tensors
    /// @param Ef Fermi energy
    /// @param T Temperature
    /// @return Tuple of {sigma, alpha, kappa}
    std::tuple<Eigen::Matrix3d, Eigen::Matrix3d, Eigen::Matrix3d> computeTransportTensors(double Ef, double T);


private:
    const Hamiltonian& H;
    const Mesh& mesh;
    double eta;
    bool temperature_in_kelvin;
    double energy_scale;


    // Reusable buffers
    mutable Eigen::VectorXd evals;
    mutable Eigen::MatrixXcd evecs;
    mutable std::array<Eigen::MatrixXcd, 3> velocity_matrices;
    mutable Eigen::MatrixXcd dH_buffer;
};
