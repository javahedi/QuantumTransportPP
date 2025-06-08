#pragma once
#include "hamiltonian.hpp"
#include <Eigen/Dense>
#include <complex>
#include <cmath>

// Reusable buffer to avoid repeated allocations
namespace BerryCurvature {
    thread_local static Eigen::MatrixXcd u0, ux, uy, uxy;
    thread_local static Eigen::VectorXd evals;
    thread_local static Eigen::MatrixXcd H_plus2, H_plus1, H_minus1, H_minus2;
}

    


/**
 * @brief Calculate the Berry curvature at a given k-point using the Fukui-Hatsugai-Suzuki discretized formula
 * @param H Hamiltonian object (must implement eigensystem)
 * @param k Wavevector (k-point) in reciprocal space
 * @param dk Small displacement in k-space
 * @param band_index Index of the band to compute curvature for
 * @return Berry curvature at the specified k-point
 */
inline double berryCurvatureFHS(const Hamiltonian& H, 
    const Eigen::Vector3d& k, double dk = 1e-3, int band_index = 0) {
        
    using namespace BerryCurvature;
    
    const Eigen::Vector3d dkx(dk, 0, 0);
    const Eigen::Vector3d dky(0, dk, 0);
    const Eigen::Vector3d dkxy(dk, dk, 0);


    //Eigen::MatrixXcd u0, ux, uy, uxy;
    //Eigen::VectorXd evals;

    // Eigenstates at the four corners of the plaquette
    H.eigensystem(k, evals, u0);
    H.eigensystem(k + dkx, evals, ux);
    H.eigensystem(k + dkxy, evals, uxy);
    H.eigensystem(k + dky, evals, uy);

    //auto overlap = [](const Eigen::VectorXcd& a, const Eigen::VectorXcd& b) {
    //    double norm = a.norm() * b.norm();
    //    return (norm < 1e-10) ? 0.0 : (a.adjoint() * b)(0);
    //};

    // Link variables (complex overlaps)
    const std::complex<double> U1 = u0.col(band_index).adjoint() * ux.col(band_index);
    const std::complex<double> U2 = ux.col(band_index).adjoint() * uxy.col(band_index);
    const std::complex<double> U3 = uxy.col(band_index).adjoint() * uy.col(band_index);
    const std::complex<double> U4 = uy.col(band_index).adjoint() * u0.col(band_index);

    
    const double product_magnitude = std::abs(U1 * U2 * U3 * U4);
    if (product_magnitude < 1e-12) return 0.0;


    double phase = std::arg(U1 * U2 * U3 * U4);
    phase = std::remainder(phase, 2 * M_PI);  // Ensure phase ∈ [-π, π]

    return phase / (dk * dk);
    
}



/// Numerically approximate ∂H/∂kx and ∂H/∂ky using finite difference
//inline Eigen::MatrixXcd dHdk(const Hamiltonian& H, const Eigen::Vector3d& k, int direction, double dk) {
//    using Vec = Eigen::Vector3d;
//    Eigen::MatrixXcd H_plus = H.Hk(k + dk * Vec(direction == 0, direction == 1, 0));
//    Eigen::MatrixXcd H_minus = H.Hk(k - dk * Vec(direction == 0, direction == 1, 0));
//    return (H_plus - H_minus) / (2.0 * dk);
//}

inline Eigen::MatrixXcd dHdk(const Hamiltonian& H, const Eigen::Vector3d& k, 
                            int direction, double dk) {

    using namespace BerryCurvature;
    
    Eigen::Vector3d dk_vec = Eigen::Vector3d::Zero();
    dk_vec(direction) = dk;
    
    // Central difference with 4th-order accuracy
    H_plus2 = H.Hk(k + 2.0 * dk_vec);
    H_plus1 = H.Hk(k + dk_vec);
    H_minus1 = H.Hk(k - dk_vec);
    H_minus2 = H.Hk(k - 2.0 * dk_vec);
    
    return (-H_plus2 + 8.0 * H_plus1 - 8.0 * H_minus1 + H_minus2) / (12.0 * dk);
}

/**
 * @brief Calculate the Berry curvature at a given k-point using the differential method
 * @param H Hamiltonian object (must implement eigensystem)
 * @param k Wavevector (k-point) in reciprocal space
 * @param dk Small displacement in k-space
 * @param band_index Index of the band to compute curvature for
 * @return Berry curvature at the specified k-point
 */
inline double berryCurvatureDifferential(const Hamiltonian& H, 
                                        const Eigen::Vector3d& k,  
                                        double dk = 1e-3, 
                                        int band_index = 0) {

    using namespace BerryCurvature;

    //Eigen::VectorXd evals;
    //Eigen::MatrixXcd evecs;

    H.eigensystem(k, evals, u0); // using u0 as a reusable buffer

    // Compute numerical derivatives ∂H/∂kx and ∂H/∂ky
    Eigen::MatrixXcd dHdkx = dHdk(H, k, 0, dk);
    Eigen::MatrixXcd dHdky = dHdk(H, k, 1, dk);

    std::complex<double> omega = 0.0;
    const auto& u_n = u0.col(band_index);
    double E_n = evals(band_index);

    for (int m = 0; m < u0.cols(); ++m) {
        if (m == band_index) continue;

        const double E_m = evals(m);
        const double deltaE = E_n - E_m;

        if (std::abs(deltaE) < 1e-8) continue;  // skip nearly degenerate bands


        const std::complex<double> vx_nm = u_n.adjoint() * dHdkx * u0.col(m);
        const std::complex<double> vy_mn = u0.col(m).adjoint() * dHdky * u_n;

        
        //double epsi = 1e-6;  // small positive number to avoid division by zero
        //double denom = std::pow(E_n - E_m, 2) + epsi;  // ε = 1e-6 or smaller
        omega += vx_nm * vy_mn / (deltaE * deltaE); // Use complex division
    }
 

    double curvature = -2.0 * std::imag(omega);

    const double BMAX = 1e2;
    if (curvature < -BMAX) curvature = -BMAX;
    if (curvature >  BMAX) curvature =  BMAX;
    
    return curvature;

    
}
