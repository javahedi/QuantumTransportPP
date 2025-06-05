#pragma once

#include "hamiltonian.hpp"
#include <complex>
#include <Eigen/Dense>
#include <cmath>

using namespace Eigen;
using std::complex;


/// @brief Altermagnet model for a two-dimensional system
/// @file altermagnet.hpp
class AltermagnetModel : public Hamiltonian {
public:
    /// @brief Constructor with default parameters
    /// @param t Hopping amplitude
    /// @param J Spin splitting parameter
    /// @param lambda Spin-orbit coupling strength
    AltermagnetModel(double t = 1.0, double J = 0.1, double lambda = 0.2)
        : t(t), J(J), lambda(lambda) {}

    /// @brief Return the Hamiltonian matrix H(k) at wavevector k
    /// @details The Hamiltonian includes scalar hopping, altermagnetic anisotropic spin splitting,
    Mat Hk(const Vec& k) const override {
        using namespace std::complex_literals;
        using std::cos, std::sin;

        // Extract components of the wavevector k
        double kx = k(0), ky = k(1);
    

        // Scalar hopping (times identity matrix)
        double eps = -2.0 * t * (cos(kx) + cos(ky));

        // Altermagnetic anisotropic spin splitting
        double d_z = J * (cos(kx) - cos(ky));

        // Optional spin-orbit terms
        double d_x = lambda * sin((kx + ky) / 2.0);
        double d_y = lambda * sin((ky - kx) / 2.0);

        // Build Hamiltonian: H = eps * σ₀ + d ⋅ σ
        Mat H = Mat::Zero(2, 2);  // Initialize a 2x2 matrix

        H += eps * Mat::Identity(2, 2); // Scalar 
        H += d_x * sigma_x();
        H += d_y * sigma_y();
        H += d_z * sigma_z();
        return H;

    }

private:
    double t;      // Hopping amplitude
    double J;      // Spin splitting parameter
    double lambda; // Spin-orbit coupling strength

    /// @brief Pauli matrices \sigma_x
    static Mat sigma_x() {
        Mat m(2, 2);
        m << 0, 1,
              1, 0;
        return m;
    }

    /// @brief Pauli matrices \sigma_y
    static Mat sigma_y() {
        Mat m(2, 2);
        m << 0, -std::complex<double>(0,1),
              std::complex<double>(0,1), 0;
        return m;
    }

    /// @brief Pauli matrices \sigma_z
    static Mat sigma_z() {
        Mat m(2, 2);
        m << 1, 0,
              0, -1;
        return m;
    }

};