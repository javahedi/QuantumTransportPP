
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>  // for std::setprecision

// ANSI escape codes for colors
const std::string RESET = "\033[0m";
const std::string BLUE = "\033[34m";
const std::string GREEN = "\033[32m";
const std::string RED = "\033[31m";

// Simple logging function
void logInfo(const std::string& msg) {
    std::cout << BLUE << "[INFO] " << RESET << msg << std::endl;
}
void logError(const std::string& msg) {
    std::cerr << RED << "[ERROR] " << RESET << msg << std::endl;
}

void logSuccess(const std::string& msg) {
    std::cout << GREEN << "[SUCCESS] " << RESET << msg << std::endl;
}