#include "Data.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>

Data::Data(const char* filename)
{
    std::fstream fin(filename, std::ios::in);

    if(!fin)
    {
        std::cerr << "Couldn't open data file " << filename << "." << std::endl;
        return;
    }

    // Get header line
    std::string line;
    std::getline(fin, line);

    while(std::getline(fin, line))
    {
        std::stringstream ss(line);
        std::string value;

        auto read_double = [&](std::stringstream& ss)
        {
            std::string temp;
            std::getline(ss, temp, ',');
            return std::stod(temp);
        };

        // then:
        x.push_back(read_double(ss));
        y.push_back(read_double(ss));
        v.push_back(read_double(ss));
        verr.push_back(read_double(ss));
        metallicity.push_back(read_double(ss));

        // Skip two values
        read_double(ss);
        read_double(ss);

        metallicity_err.push_back(read_double(ss));

//        std::cout << x.back() << ' ' << y.back() <<' ' << v.back() << ' ';
//        std::cout << verr.back() << ' ' << metallicity.back() << std::endl;
    }

    // Compute theta
    theta.resize(x.size());
    for(size_t i=0; i<x.size(); ++i)
        theta[i] = atan2(y[i], x[i]);

    std::cout << "# Loaded " << x.size() << " data points from file ";
    std::cout << filename << "." << std::endl;

    fin.close();
}

