#include "hamiltonian.hpp"
#include <stdexcept>

/// @brief  Hamiltonian class implementation
/// @file hamiltonian.cpp
/// This file implements the eigensystem method for the Hamiltonian class.
/// The eigensystem method computes the eigenvalues and eigenvectors of the Hamiltonian matrix H(k)
/// at a given wavevector k using Eigen's SelfAdjointEigenSolver.
/// @details The Hamiltonian matrix is expected to be Hermitian, and the method will throw an exception
/// if the eigensystem computation fails.
/// @param k 
/// @param evals 
/// @param evecs 
void Hamiltonian::eigensystem(const Vec& k, Eigen::VectorXd& evals, Eigen::MatrixXcd& evecs) const {
    // Compute the Hamiltonian matrix at wavevector k
    Mat Hk_matrix = Hk(k);
    
    // Use Eigen's SelfAdjointEigenSolver for Hermitian matrices
    Eigen::SelfAdjointEigenSolver<Mat> solver(Hk_matrix);
    
    if (solver.info() != Eigen::Success) {
        // Handle the case where the eigensystem computation fails
        throw std::runtime_error("Eigensystem computation failed");
    }
    
    evals = solver.eigenvalues();
    evecs = solver.eigenvectors();
}