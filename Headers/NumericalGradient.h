#pragma once
#ifndef NUMERICAL_GRADIENT_H
#define NUMERICAL_GRADIENT_H

#include <functional>
#include <vector>

class NumericalGradient
{
public:
    using Evaluator = std::function<double(const std::vector<double> &)>;

    NumericalGradient(Evaluator evaluator, double baseStep);

    std::vector<double> gradient(const std::vector<double> &z) const;

private:
    Evaluator F;
    double h;
};

#endif
