#pragma once
#ifndef POLYNOMIAL_OBJECTIVE_H
#define POLYNOMIAL_OBJECTIVE_H

#include "ObjectiveFunction.h"
#include <vector>

class PolynomialObjective : public ObjectiveFunction
{
public:
    PolynomialObjective(
        const std::vector<double>& x,
        const std::vector<double>& y);

    double evaluate(
        const std::vector<double>& theta) const override;
    
    double polynomial(
        const std::vector<double>& theta,
        double x) const;

    void printPolynomial(const std::vector<double>& theta) const;
    
private:
    std::vector<double> xData;
    std::vector<double> yData;
};

// Degree-by-degree polynomial fitting using the project's gradient descent
// optimizer. Returns coefficients in ascending powers of x.
std::vector<double> Polynomial_Approximator_Gradient_Descent(
    const std::vector<double>& x,
    const std::vector<double>& y);

#endif
