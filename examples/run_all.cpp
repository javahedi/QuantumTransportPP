// src/main.cpp or examples/run_all.cpp
#include "haldane.hpp"
#include "mesh.hpp"
#include "dos.hpp"
#include <iostream>

int main() {
    Mesh mesh(30, 30);
    HaldaneModel model(1.0, 0.1, M_PI/2);

    DOS dos(model, mesh);
    auto spectrum = dos.computeDOS(-1.0, 1.0, 50, 1e-2);

    std::cout << "DOS computed for Haldane model.\n";
    return 0;
}
