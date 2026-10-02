#include "PolynomialObjective.h"

#include <fstream>
#include <iostream>
#include <stdexcept>

int main()
{
    try
    {
        std::ifstream input("dataset.txt");
        if (!input)
            throw std::runtime_error("Unable to open dataset.txt.");

        std::vector<double> x, y;
        double xi, yi;
        while (input >> xi >> yi)
        {
            x.push_back(xi);
            y.push_back(yi);
        }
        if (x.empty())
            throw std::runtime_error("No data found in dataset.txt.");

        Polynomial_Approximator_Gradient_Descent(x, y);
    }
    catch (const std::exception& error)
    {
        std::cerr << "Gradient descent polynomial approximation failed: "
                  << error.what() << '\n';
        return 1;
    }
    return 0;
}
