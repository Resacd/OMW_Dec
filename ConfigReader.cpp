#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include "ConfigReader.h"


/*
ConfigReader::ConfigReader(const std::string& configFile) : configFile(configFile) {}


bool ConfigReader::readConfigData(std::map<std::string, double>& configData) {
    std::ifstream file(configFile);

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << configFile << std::endl;
        return false;
    }

    std::string line;
    while (std::getline(file, line)) {
        // Skip comment lines and empty lines
        if (line.empty() || line[0] == '#') {
            continue;
        }

        size_t pos = line.find_first_of(" \t");
        if (pos != std::string::npos) {
            std::string key = line.substr(0, pos);
            std::string value = line.substr(pos + 1);

            try {
                configData[key] = std::stod(value);
            }
            catch (const std::invalid_argument& e) {
                std::cerr << "Error converting value to double for key " << key << ": " << e.what() << std::endl;
                file.close();
                return false;
            }
            catch (const std::out_of_range& e) {
                std::cerr << "Error converting value to double (out of range) for key " << key << ": " << e.what() << std::endl;
                file.close();
                return false;
            }
        }
    }

    file.close();
    return true;
}


bool ConfigReader::readEquations(std::vector<std::pair<std::string, std::string>>& equations) {
    std::ifstream file(configFile);

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << configFile << std::endl;
        return false;
    }

    std::string line;
    std::string equation;
    std::string conditions;

    while (std::getline(file, line)) {
        // Skip comment lines and empty lines
        if (line.empty() || line[0] == '#') {
            continue;
        }

        size_t equationPos = line.find("Equation");
        size_t conditionsPos = line.find("Conditions:");

        if (equationPos != std::string::npos) {
            if (!equation.empty() && !conditions.empty()) {
                equations.emplace_back(conditions, equation);
                equation.clear();
                conditions.clear();
            }
        }
        else if (conditionsPos != std::string::npos) {
            conditions = line.substr(conditionsPos + strlen("Conditions:"));
        }
        else {
            equation += line + '\n';
        }
    }

    if (!equation.empty() && !conditions.empty()) {
        equations.emplace_back(conditions, equation);
    }

    file.close();
    return true;
}


void ConfigReader::printEquations(const std::vector<std::pair<std::string, std::string>>& equations) const {
    for (const auto& eq : equations) {
        std::cout << "Equation Section: " << eq.first << "\n" << eq.second << "\n\n";
    }
}
*/
