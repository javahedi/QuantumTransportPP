#pragma once
#include "hamiltonian.hpp"
#include <Eigen/Dense>
#include <complex>
#include <cmath>

/*
(Fukui-Hatsugai-Suzuki) or discretized Berry curvature.
*/
inline double berryCurvature(const Hamiltonian& H, const Eigen::Vector3d& k, double dk = 1e-3, int band_index = 0) {
    using Vec = Eigen::Vector3d;
    using std::arg, std::abs;

    Eigen::MatrixXcd u0, ux, uy, uxy;
    Eigen::VectorXd evals;

    // Eigenstates at the four corners of the plaquette
    H.eigensystem(k, evals, u0);
    H.eigensystem(k + Vec(dk, 0, 0), evals, ux);
    H.eigensystem(k + Vec(dk, dk, 0), evals, uxy);
    H.eigensystem(k + Vec(0, dk, 0), evals, uy);

    // Link variables (U(1) Wilson lines)
    auto U1 = (u0.col(band_index).adjoint() * ux.col(band_index)).normalized();
    auto U2 = (ux.col(band_index).adjoint() * uxy.col(band_index)).normalized();
    auto U3 = (uxy.col(band_index).adjoint() * uy.col(band_index)).normalized();
    auto U4 = (uy.col(band_index).adjoint() * u0.col(band_index)).normalized();

    // Wilson loop and Berry curvature
    std::complex<double> loop = U1 * U2 * U3 * U4;
    double curvature = arg(loop) / (dk * dk);

    return curvature;
}
