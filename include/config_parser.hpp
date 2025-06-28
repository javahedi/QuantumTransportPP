#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

class ConfigParser {
public:
    explicit ConfigParser(const std::string& filename) {
        parseFile(filename);
    }

    std::vector<double> getRange(const std::string& section, 
                               const std::string& key,
                               const std::vector<double>& defaults = {}) const {
        try {
            std::string val = m_sections.at(section).at(key);
            if (val.find(':') != std::string::npos) {
                return parseRange(val);
            }
            return parseList(val);
        } catch (...) {
            if (defaults.empty()) throw;
            return defaults;
        }
    }

    template<typename T>
    T get(const std::string& section, 
          const std::string& key,
          T default_value) const {
        try {
            std::string val = m_sections.at(section).at(key);
            std::istringstream iss(val);
            T result;
            iss >> result;
            return result;
        } catch (...) {
            return default_value;
        }
    }

private:
    std::unordered_map<std::string, 
                      std::unordered_map<std::string, std::string>> m_sections;

    void parseFile(const std::string& filename) {
        std::ifstream file(filename);
        if (!file) throw std::runtime_error("Config file not found");

        std::string current_section;
        std::string line;

        while (std::getline(file, line)) {
            line = trim(line);
            if (line.empty() || line[0] == ';' || line[0] == '#') continue;

            if (line[0] == '[' && line.back() == ']') {
                current_section = line.substr(1, line.size() - 2);
            } else {
                size_t eq_pos = line.find('=');
                if (eq_pos != std::string::npos) {
                    std::string key = trim(line.substr(0, eq_pos));
                    std::string value = trim(line.substr(eq_pos + 1));
                    m_sections[current_section][key] = value;
                }
            }
        }
    }

    static std::vector<double> parseRange(const std::string& str) {
        std::istringstream iss(str);
        double start, end, step;
        char colon;
        if (!(iss >> start >> colon >> end >> colon >> step) || colon != ':') {
            throw std::runtime_error("Invalid range format");
        }
        
        std::vector<double> result;
        for (double val = start; val <= end + 1e-9; val += step) {
            result.push_back(val);
        }
        return result;
    }

    static std::vector<double> parseList(const std::string& str) {
        std::vector<double> result;
        std::istringstream iss(str);
        std::string token;
        while (std::getline(iss, token, ',')) {
            result.push_back(std::stod(trim(token)));
        }
        return result;
    }

    static std::string trim(const std::string& s) {
        size_t start = s.find_first_not_of(" \t");
        size_t end = s.find_last_not_of(" \t");
        return (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
    }
};