
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

    ///  @brief Compute real part of conductivity tensor σ_ij
    /// @param Ef Fermi energy
    /// @param T Temperature in Kelvin (if temperature_in_kelvin is true)
    /// @return Conductivity tensor σ_ij as a 3x3 matrix
    Eigen::Matrix3d conductivity(double Ef, double T);

private:
    const Hamiltonian& H;
    const Mesh& mesh;
    double eta;
    bool temperature_in_kelvin;
    double energy_scale;
};
