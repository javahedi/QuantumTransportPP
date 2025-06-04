cmake_minimum_required(VERSION 3.14)
project(QuantumTransport++ LANGUAGES CXX)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


# set the include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Source files
file(GLOB_RECURSE SOURCES src/*.cpp)


# Create a shared library from your source files (if reusable)
add_library(QtTransportLib STATIC ${SOURCES})


# Example executables for demos
add_executable(run_haldane examples/run_haldane.cpp)
target_link_libraries(run_haldane QtTransportLib)

# Test executables
add_executable(test_mesh tests/test_mesh.cpp)
target_link_libraries(test_mesh QtTransportLib)

add_executable(test_dos tests/test_dos.cpp)
target_link_libraries(test_dos QtTransportLib)

add_executable(test_bandstructure tests/test_bandstructure.cpp)
target_link_libraries(test_bandstructure QtTransportLib)

add_executable(test_berry tests/test_berry.cpp)
target_link_libraries(test_berry QtTransportLib)

add_executable(test_boltzmann tests/test_boltzmann.cpp)
target_link_libraries(test_boltzmann QtTransportLib)

add_executable(test_kubo tests/test_kubo.cpp)
target_link_libraries(test_kubo QtTransportLib)

# ===== Test Executables =====
#add_executable(test_mesh tests/test_mesh.cpp ${SOURCES})
#add_executable(test_bandstructure tests/test_bandstructure.cpp ${SOURCES})
#add_executable(test_berry tests/test_berry.cpp ${SOURCES})
#add_executable(test_boltzmann tests/test_boltzmann.cpp ${SOURCES})
#add_executable(test_kubo tests/test_kubo.cpp ${SOURCES})
#add_executable(test_dos tests/test_dos.cpp ${SOURCES})

