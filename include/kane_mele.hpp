#pragma once
#include "hamiltonian.hpp"
#include <complex>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>

class KaneMeleModel : public Hamiltonian {
public:

    double t = 1.0;         // Hopping parameter
    double lambda_SO = 0.1; // Spin-orbit coupling strength
    double lambda_v = 0.2 ; // sublattice potential
    bool include_Rashba = false;
    double lambda_R = 0.0;

    KaneMeleModel(double t_=1.0, double so_=0.1, double v_=0.2, bool rashba=false) 
        : t(t_), lambda_SO(so_), lambda_v(v_), include_Rashba(rashba) {}

    
    Mat Hk(const Vec& k) const override {
        using namespace std::complex_literals; 
        double kx = k(0), ky = k(1);
        const std::complex<double> I(0, 1);

        // Pauli matrices
        Eigen::Matrix2cd sigma_x, sigma_y, sigma_z;
        sigma_x << 0, 1,
                   1, 0;
        sigma_y << 0, -I,
                   I, 0;
        sigma_z << 1, 0,
                   0, -1;
                   
        // Sublattice part (similar to graphene)
        //double f_re =  std::cos(kx) + 2.0 * std::cos(kx / 2.0) * std::cos(std::sqrt(3) * ky / 2.0);
        //double f_im = -std::sin(kx) - 2.0 * std::sin(kx / 2.0) * std::cos(std::sqrt(3) * ky / 2.0);
        //std::complex<double> f = f_re + I * f_im;
        std::complex<double> f = std::exp(I*kx) 
                       + 2.0 * std::exp(I*kx/2.0) * std::cos(std::sqrt(3)*ky/2.0);

        // Build 2x2 H(k) ignoring spin
        Eigen::Matrix2cd H0;
        H0 = -t * f.real() * sigma_x - t * f.imag() * sigma_y + lambda_v * sigma_z;


        // Spin-orbit term: 2x2 σ_z ⊗ s_z => expands to 4x4
        Eigen::Matrix4cd H = Eigen::Matrix4cd::Zero();
        H.block<2,2>(0,0) = H0 + lambda_SO * sigma_z;  // spin up
        H.block<2,2>(2,2) = H0 - lambda_SO * sigma_z;   // spin down

        // Verify time-reversal symmetry: TH(k)T⁻¹ = H(-k)
        Eigen::Matrix4cd T = Eigen::kroneckerProduct(sigma_x, Eigen::Matrix2cd::Identity());
        T *= I;  // T = iσ_x ⊗ I
        assert((T * Hk(k) * T.inverse() - Hk(-k)).norm() < 1e-10);

        // Rashba coupling (if needed)
        if (include_Rashba) {
            H.block<2,2>(0,2) = lambda_R*(std::sin(kx)*sigma_x - std::sin(ky)*sigma_y);
            H.block<2,2>(2,0) = H.block<2,2>(0,2).adjoint();
        }
        

        return H;
    }

};