#pragma once
#include <Eigen/Dense>


// Abstract base class for band structure calculations
/// @file bandstructure.hpp
/// @brief Abstract base class for band structure calculations
class BandStructure {
public:
    using Vec = Eigen::Vector3d;

    /// @brief Virtual destructor for proper cleanup of derived classes
    virtual ~BandStructure() = default;

    /// @brief Energy dispersion ε(k)
    virtual double energy(const Vec& k) const = 0;


    /// @brief Group velocity ∇k ε(k)
    virtual Vec groupVelocity(const Vec& k) const = 0;
    
};