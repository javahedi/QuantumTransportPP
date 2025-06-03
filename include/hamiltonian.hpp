#pragma once
#include <Eigen/Dense>

// Abstract base class for Hamiltonian models
class Hamiltonian {
public:
    using Vec = Eigen::Vector3d;
    using Mat = Eigen::MatrixXcd;

    virtual ~Hamiltonian() = default;

    virtual Mat Hk(const Vec& k) const = 0; // Return complex matrix H(k)

    virtual void eigensystem(const Vec& k, Eigen::VectorXd& evals, Eigen::MatrixXcd& evecs) const; 
    
};