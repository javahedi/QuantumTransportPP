#include "mesh.hpp"
#include <iostream>


int main() {
    
    Mesh mesh(5,5); // 2D 5x5 mesh
    const auto& kpoints = mesh.getKPoints();

    std::cout << "Generated " << mesh.size() << " k-points:\n";
    for (const auto& kp : kpoints) {
        std::cout << kp.transpose() << "\n";
    }

    return 0;
}

/*
mkdir build && cd build
cmake ..
make
./test_mesh
*/