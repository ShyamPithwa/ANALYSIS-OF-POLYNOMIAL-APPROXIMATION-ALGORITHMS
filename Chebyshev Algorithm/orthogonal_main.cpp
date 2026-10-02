#include "factorizations.h"

#include <iostream>
#include <stdexcept>

int main()
{
    try
    {
        std::vector<double> x, y;
        readDataPoints("dataset.txt", x, y);
        Orthogonal_Polynomial_Approximator(x, y);
    }
    catch (const std::exception& error)
    {
        std::cerr << "Orthogonal polynomial approximation failed: "
                  << error.what() << '\n';
        return 1;
    }
    return 0;
}
