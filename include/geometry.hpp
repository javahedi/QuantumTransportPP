#pragma once
#include "hamiltonian.hpp"
#include <Eigen/Dense>
#include <complex>
#include <cmath>

/**
 * @brief Calculate the Berry curvature at a given k-point using the Fukui-Hatsugai-Suzuki discretized formula
 * @param H Hamiltonian object (must implement eigensystem)
 * @param k Wavevector (k-point) in reciprocal space
 * @param dk Small displacement in k-space
 * @param band_index Index of the band to compute curvature for
 * @return Berry curvature at the specified k-point
 */
inline double berryCurvatureFHS(const Hamiltonian& H, const Eigen::Vector3d& k, double dk = 1e-3, int band_index = 0) {
    using Vec = Eigen::Vector3d;
    using std::arg;

    Eigen::MatrixXcd u0, ux, uy, uxy;
    Eigen::VectorXd evals;

    // Eigenstates at the four corners of the plaquette
    H.eigensystem(k, evals, u0);
    H.eigensystem(k + Vec(dk, 0, 0), evals, ux);
    H.eigensystem(k + Vec(dk, dk, 0), evals, uxy);
    H.eigensystem(k + Vec(0, dk, 0), evals, uy);

    // Link variables (complex overlaps)
    std::complex<double> U1 = u0.col(band_index).adjoint() * ux.col(band_index);
    std::complex<double> U2 = ux.col(band_index).adjoint() * uxy.col(band_index);
    std::complex<double> U3 = uxy.col(band_index).adjoint() * uy.col(band_index);
    std::complex<double> U4 = uy.col(band_index).adjoint() * u0.col(band_index);

    // 🔒 Safety: avoid unreliable overlaps (near orthogonal eigenstates)
    double eps = 1e-12;
    if (std::abs(U1) < eps || std::abs(U2) < eps || std::abs(U3) < eps || std::abs(U4) < eps)
        return 0.0;

    // 🔒 Optional: skip problematic symmetry points
    //if ((k.array().abs() > M_PI - 1e-6).any()) return 0.0;


    // Wilson loop and Berry curvature
    std::complex<double> loop = U1 * U2 * U3 * U4;
    double phase = arg(loop);

    // Optional unwrap (in practice arg() is already in [-π, π])
    if (phase > M_PI) phase -= 2 * M_PI;
    if (phase < -M_PI) phase += 2 * M_PI;

    return phase / (dk * dk);

    double curvature = phase / (dk * dk);
    const double BMAX = 1e2;
    if (std::abs(curvature) > BMAX)
        curvature = (curvature > 0 ? BMAX : -BMAX);
    return curvature;
    
}



/// Numerically approximate ∂H/∂kx and ∂H/∂ky using finite difference
inline Eigen::MatrixXcd dHdk(const Hamiltonian& H, const Eigen::Vector3d& k, int direction, double dk) {
    using Vec = Eigen::Vector3d;
    Eigen::MatrixXcd H_plus = H.Hk(k + dk * Vec(direction == 0, direction == 1, 0));
    Eigen::MatrixXcd H_minus = H.Hk(k - dk * Vec(direction == 0, direction == 1, 0));
    return (H_plus - H_minus) / (2.0 * dk);
}



inline double berryCurvatureDifferential(const Hamiltonian& H, const Eigen::Vector3d& k,  double dk = 1e-3, int band_index = 0) {
    using Mat = Eigen::MatrixXcd;
    using Vec = Eigen::Vector3d;

    Eigen::MatrixXcd Hk = H.Hk(k);
    Eigen::VectorXd evals;
    Eigen::MatrixXcd evecs;

    H.eigensystem(k, evals, evecs);

    // Compute numerical derivatives ∂H/∂kx and ∂H/∂ky
    Mat dHdkx = dHdk(H, k, 0, dk);
    Mat dHdky = dHdk(H, k, 1, dk);

    std::complex<double> omega = 0.0;
    const auto& u_n = evecs.col(band_index);
    double E_n = evals(band_index);

    for (int m = 0; m < evecs.cols(); ++m) {
        if (m == band_index) continue;
        const auto& u_m = evecs.col(m);
        double E_m = evals(m);



        std::complex<double> vx_nm = u_n.adjoint() * dHdkx * u_m;
        std::complex<double> vy_mn = u_m.adjoint() * dHdky * u_n;

        //if (std::abs(E_n - E_m) < 1e-6) continue;  // skip nearly degenerate bands
        double epsi = 1e-6;  // small positive number to avoid division by zero
        double denom = std::pow(E_n - E_m, 2) + epsi;  // ε = 1e-6 or smaller
        omega += vx_nm * vy_mn / std::pow(E_n - E_m, 2);
    }
 

    double curvature = -2.0 * std::imag(omega);
    const double BMAX = 1e2;
    if (std::abs(curvature) > BMAX)
        curvature = (curvature > 0 ? BMAX : -BMAX);
    return curvature;

    
}
