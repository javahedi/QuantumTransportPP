#pragma once
#include <Eigen/Dense>


// Abstract base class for band structure calculations
class BandStructure {
public:
    using Vec = Eigen::Vector3d;

    // Virtual destructor for proper cleanup of derived classes
    virtual ~BandStructure() = default;

    // Energy dispersion ε(k)
    virtual double energy(const Vec& k) const = 0;


    // Group velocity ∇k ε(k)
    virtual Vec groupVelocity(const Vec& k) const = 0;
    
};