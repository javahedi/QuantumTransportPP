cmake_minimum_required(VERSION 3.14)
project(QuantumTransport++ LANGUAGES CXX)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# set the include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# set the source files
set(SOURCES
    src/mesh.cpp
    src/hamiltonian.cpp
    src/boltzmann.cpp
    src/kubo.cpp
)

# add test executable
add_executable(test_mesh tests/test_mesh.cpp ${SOURCES})
add_executable(test_bandstructure tests/test_bandstructure.cpp ${SOURCES})
add_executable(test_berry tests/test_berry.cpp ${SOURCES})
add_executable(test_boltzmann tests/test_boltzmann.cpp ${SOURCES})
add_executable(test_kubo tests/test_kubo.cpp ${SOURCES})

