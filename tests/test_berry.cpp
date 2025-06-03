#include "mesh.hpp"
#include "haldane.hpp"
#include "geometry.hpp"
#include <iostream>


int main() {
    HaldaneModel model;
    Mesh mesh(10, 10);

    for (const auto& k : mesh.getKPoints()) {
        double Omega = berryCurvature(model, k);
        std::cout << k.head<2>().transpose() << " " << Omega << std::endl;
    }
}
