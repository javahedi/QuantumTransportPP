#pragma once
#include "hamiltonian.hpp"
#include "mesh.hpp"
#include "geometry.hpp"
#include <Eigen/Dense>

// Solver for semiclassical Boltzmann transport with optional quantum geometry effects
class BoltzmannSolver {
public:
    // Constructor
    // @param H: Hamiltonian model (must implement eigensystem)
    // @param mesh: k-point mesh in BZ
    // @param tau: relaxation time
    // @param temperature_in_kelvin: if true, T is interpreted in Kelvin
    // @param energy_scale: energy unit (e.g., t1), used when converting k_B T
    BoltzmannSolver(const Hamiltonian& H, const Mesh& mesh, double tau,
                    bool temperature_in_kelvin = false,
                    double energy_scale = 1.0)
        : H(H), mesh(mesh), tau(tau),
          temperature_in_kelvin(temperature_in_kelvin),
          energy_scale(energy_scale) {}

    // Compute group velocity and anomalous velocity: v + Omega x E
    Eigen::Vector3d velocity(const Eigen::Vector3d& k, int band, const Eigen::Vector3d& E, double dk = 1e-4) const;

    // Compute conductivity tensor σ_ij via Boltzmann (RTA)
    Eigen::Matrix3d conductivity(double Ef, double T, const Eigen::Vector3d& Efield);

private:
    const Hamiltonian& H;
    const Mesh& mesh;
    double tau;
    bool temperature_in_kelvin;
    double energy_scale;
};
