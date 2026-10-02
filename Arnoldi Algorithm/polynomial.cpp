#include "Arnoldi_iteration.h"
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H
#include<iostream>

#include <cmath>
#include <algorithm>

Polynomial::Polynomial() {
    coeff.push_back(0.0);
}

Polynomial::Polynomial(const vector<double>& coefficients) {
    coeff = coefficients;
}

// p(x) -> x * p(x)
void Polynomial::increasePower() {

    coeff.insert(coeff.begin(), 0.0);
}

// p(x) -> p(x) + alpha * q(x)
void Polynomial::addPolynomial(
    const Polynomial& q,
    double alpha
) {

    const std::size_t maxSize = std::max(coeff.size(), q.coeff.size());

    coeff.resize(maxSize, 0.0);

    for (std::size_t i = 0; i < q.coeff.size(); ++i) {
        coeff[i] += alpha * q.coeff[i];
    }
}

// p(x) -> alpha * p(x)
void Polynomial::scale(double alpha) {

    for (double& c : coeff) {
        c *= alpha;
    }
}

// Evaluate p(x)
double Polynomial::evaluate(double x) const {

    double result = 0.0;

    for (std::size_t i = 0; i < coeff.size(); ++i) {
        result += coeff[i] * pow(x, i);
    }

    return result;
}

// Return coefficients
vector<double> Polynomial::getCoefficients() const {
    return coeff;
}

// Return degree
int Polynomial::degree() const {
    return coeff.size() - 1;
}

#endif