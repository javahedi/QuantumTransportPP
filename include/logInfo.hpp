// include/logInfo.hpp
#ifndef LOG_INFO_HPP
#define LOG_INFO_HPP

#include <iostream>
#include <string>
#include <chrono>
#include <iomanip>

/** \file logInfo.hpp
 * \brief Provides colored logging functions for console output.
 * \details Defines logInfo, logError, and logSuccess for status messages with ANSI color codes.
 * \author [Your Name]
 * \date June 2025
 */

namespace QuantumTransport {
    // ANSI escape codes for colors
    const std::string RESET = "\033[0m";
    const std::string BLUE = "\033[34m";
    const std::string GREEN = "\033[32m";
    const std::string RED = "\033[31m";

    /** \brief Logs an informational message to stdout in blue.
     * \param msg The message to log.
     */
    void logInfo(const std::string& msg) {
        std::cout << BLUE << "[INFO] " << RESET << msg << std::endl;
    }

    /** \brief Logs an error message to stderr in red.
     * \param msg The error message to log.
     */
    void logError(const std::string& msg) {
        std::cerr << RED << "[ERROR] " << RESET << msg << std::endl;
    }

    /** \brief Logs a success message to stdout in green.
     * \param msg The success message to log.
     */
    void logSuccess(const std::string& msg) {
        std::cout << GREEN << "[SUCCESS] " << RESET << msg << std::endl;
    }
} // namespace QuantumTransport



/// @brief class for displaying progress in console applications. 

class ProgressBar {
private:
    size_t total;
    size_t current;
    double start_time;
    bool parallel;
    
public:
    ProgressBar(size_t total, bool parallel = false) 
        : total(total), current(0), parallel(parallel) {
        start_time = std::chrono::duration<double>(
            std::chrono::steady_clock::now().time_since_epoch()).count();
    }
    
    void update(size_t n = 1) {
        if (parallel) {
            #pragma omp atomic
            current += n;
        } else {
            current += n;
        }
        
        print();
    }
    
    void print() {
        const int barWidth = 50;
        float progress = static_cast<float>(current) / total;
        double now = std::chrono::duration<double>(
            std::chrono::steady_clock::now().time_since_epoch()).count();
        double elapsed = now - start_time;
        
        std::cout << "[";
        int pos = barWidth * progress;
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << int(progress * 100.0) << "% "
                  << "(" << current << "/" << total << ") "
                  << std::fixed << std::setprecision(1) << elapsed << "s\r";
        std::cout.flush();
    }
    
    void finish() {
        current = total;
        print();
        std::cout << std::endl;
    }
};

#endif // LOG_INFO_HPP