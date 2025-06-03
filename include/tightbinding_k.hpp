#pragma once
#include "bandstructure.hpp"
#include <iostream>
#include <cmath>


class TightBindingSquare : public BandStructure {
public:
    double t;

    TightBindingSquare(double t_in = 1.0) : t(t_in) {}

    double energy(const Vec& k) const override {
        return -2.0 * t * (std::cos(k(0)) + std::cos(k(1)));
    }

    Vec groupVelocity(const Vec& k) const override {
        Vec v;
        v(0) = 2 * t * std::sin(k(0));
        v(1) = 2 * t * std::sin(k(1));
        v(2) = 0.0; // No z-component in 2D
        return v;
    }
};
