#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <cmath>
class ConfigReader {
public:

    ConfigReader(const std::string& configFile);
    bool readConfigData(std::map<std::string, double>& constants);
    bool readEquations(std::vector<std::pair<std::string, std::string>>& equations);
    void printEquations(const std::vector<std::pair<std::string, std::string>>& equations) const;

private:
    std::string configFile;
};

