cmake_minimum_required(VERSION 3.14)
project(QuantumTransport++ LANGUAGES CXX)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Source files
file(GLOB_RECURSE SOURCES src/*.cpp)

# OpenMP configuration for AppleClang with Homebrew's libomp
if(APPLE)
    # Find Homebrew's libomp
    find_path(LIBOMP_INCLUDE_DIR omp.h PATHS /opt/homebrew/opt/libomp/include REQUIRED)
    find_library(LIBOMP_LIBRARY omp PATHS /opt/homebrew/opt/libomp/lib REQUIRED)
    
    message(STATUS "Found OpenMP include dir: ${LIBOMP_INCLUDE_DIR}")
    message(STATUS "Found OpenMP library: ${LIBOMP_LIBRARY}")
    
    # Set OpenMP flags separately
    set(OPENMP_COMPILE_OPTIONS "-Xpreprocessor" "-fopenmp")
    set(OPENMP_INCLUDE_DIRS "${LIBOMP_INCLUDE_DIR}")
    set(OPENMP_LINK_LIBRARIES "${LIBOMP_LIBRARY}")
    
    include_directories(${OPENMP_INCLUDE_DIRS})
else()
    # Standard OpenMP detection for other platforms
    find_package(OpenMP REQUIRED)
    if(OpenMP_CXX_FOUND)
        set(OPENMP_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
        set(OPENMP_LINK_LIBRARIES OpenMP::OpenMP_CXX)
    endif()
endif()

# Create the library
add_library(QtTransportLib STATIC ${SOURCES})
if(APPLE OR OpenMP_CXX_FOUND)
    target_compile_options(QtTransportLib PRIVATE ${OPENMP_COMPILE_OPTIONS})
    target_link_libraries(QtTransportLib PRIVATE ${OPENMP_LINK_LIBRARIES})
endif()

# Helper function for executables
function(add_qt_executable name source)
    add_executable(${name} ${source})
    target_link_libraries(${name} PRIVATE QtTransportLib)
    if(APPLE OR OpenMP_CXX_FOUND)
        target_compile_options(${name} PRIVATE ${OPENMP_COMPILE_OPTIONS})
        target_link_libraries(${name} PRIVATE ${OPENMP_LINK_LIBRARIES})
    endif()
endfunction()


# Example executables
add_qt_executable(run_all examples/run_all.cpp)
add_qt_executable(transport_phi_haldane examples/transport_phi_haldane.cpp)
add_qt_executable(altermagnet_berry examples/altermagnet_berry.cpp)
add_qt_executable(altermagnet_conductivity examples/altermagnet_conductivity.cpp)

# Test executables
add_qt_executable(test_mesh tests/test_mesh.cpp)
add_qt_executable(test_dos tests/test_dos.cpp)
add_qt_executable(test_bandstructure tests/test_bandstructure.cpp)
add_qt_executable(test_berry tests/test_berry.cpp)
add_qt_executable(test_boltzmann tests/test_boltzmann.cpp)
add_qt_executable(test_kubo tests/test_kubo.cpp)