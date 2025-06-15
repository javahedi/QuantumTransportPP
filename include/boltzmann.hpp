#pragma once
#include "hamiltonian.hpp"
#include "mesh.hpp"
#include "geometry.hpp"
#include <Eigen/Dense>



/// @brief Struct to store velocity result and phase-space factor
struct VelocityResult {
    Eigen::Vector3d velocity;
    double phaseSpaceFactor;
};



///  @file boltzmann.hpp
///  @brief BoltzmannSolver class for calculating transport properties with with optional quantum geometry effects
class BoltzmannSolver {
public:
    // Constructor
    /// @param H: Hamiltonian model (must implement eigensystem)
    /// @param mesh: k-point mesh in BZ
    /// @param tau: relaxation time
    /// @param temperature_in_kelvin: if true, T is interpreted in Kelvin
    /// @param energy_scale: energy unit (e.g., t1), used when converting k_B T
    BoltzmannSolver(const Hamiltonian& H, const Mesh& mesh, double tau,
                    bool temperature_in_kelvin = false,
                    double energy_scale = 1.0)
        : H(H), mesh(mesh), tau(tau),
          temperature_in_kelvin(temperature_in_kelvin),
          energy_scale(energy_scale) {}

   
    /// @brief Compute velocity at k-point, including anomalous and Lorentz contributions
    /// @param k: k-point in reciprocal space
    /// @param band: band index
    /// @param E: Electric field vector
    /// @param B: Magnetic field vector
    /// @param dk: small displacement in k-space
    /// @return VelocityResult: full velocity and phase-space factor
    VelocityResult velocity(double energy, double Ef, double T,
                            const Eigen::Vector3d& k, int band, 
                            const Eigen::Vector3d& gradT,
                            const Eigen::Vector3d& E, 
                            const Eigen::Vector3d& B, 
                            double dk = 1e-4) const;

    
    double phaseSpaceFactor(const Eigen::Vector3d& k, int band, 
                            const Eigen::Vector3d& B) const;
    
    // Compute both σ_ij and α_ij in one pass
    /// @brief Compute transport tensors for conductivity and thermopower
    /// @param Ef: Fermi energy
    /// @param T: Temperature (in Kelvin if temperature_in_kelvin is true)
    /// @param Efield: Electric field vector
    /// @return tuble of conductivity tensor and thermopower tensor

    std::tuple<Eigen::Matrix3d, Eigen::Matrix3d> 
    computeTransportTensors(double Ef, double T, 
                            const Eigen::Vector3d& gradT,
                            const Eigen::Vector3d& Efield,
                            const Eigen::Vector3d& Bfield) const;


    Eigen::Matrix3d secondDerivatives(const Eigen::Vector3d& k, int band, 
                                      double dk) const;


private:
    const Hamiltonian& H;
    const Mesh& mesh;
    double tau;
    bool temperature_in_kelvin;
    double energy_scale;

    // Reusable buffers to avoid allocations
    mutable Eigen::VectorXd evals_plus;
    mutable Eigen::VectorXd evals_minus;
    mutable Eigen::MatrixXcd evecs_dummy;
    mutable Eigen::VectorXd evals;
    mutable Eigen::MatrixXcd evecs;
};
