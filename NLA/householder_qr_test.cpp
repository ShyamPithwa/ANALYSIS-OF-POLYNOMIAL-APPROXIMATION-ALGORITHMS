#include "factorizations.h"

#include <exception>
#include <iostream>
#include <vector>

int main()
{
    try
    {
        std::vector<double> x;
        std::vector<double> y;
        readDataPoints("dataset.txt", x, y);

        // Prints the degree-by-degree residuals and the best polynomial.
        Polynomial_Approximator(x, y);
    }
    catch (const std::exception& error)
    {
        std::cerr << "Householder QR polynomial approximation failed: "
                  << error.what() << '\n';
        return 1;
    }
    return 0;
}
