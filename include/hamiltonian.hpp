#pragma once
#include <Eigen/Dense>

///  @brief Hamiltonian.hpp
///  @file hamiltonian.hpp
///  @brief Abstract base class for Hamiltonian models in solid-state physics
class Hamiltonian {
public:
    using Vec = Eigen::Vector3d;
    using Mat = Eigen::MatrixXcd;

    virtual ~Hamiltonian() = default;

    /// @brief Return complex matrix H(k)
    virtual Mat Hk(const Vec& k) const = 0; // 

    /// @brief Return the eigenvalues and eigenvectors of H(k)
    virtual void eigensystem(const Vec& k, Eigen::VectorXd& evals, Eigen::MatrixXcd& evecs) const; 
    
};