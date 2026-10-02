#include "PolynomialObjective.h"
#include "GradientDescent.h"
#include "GradientDescentOptions.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>

PolynomialObjective::PolynomialObjective(
    const std::vector<double>& x,
    const std::vector<double>& y)
    : xData(x), yData(y)
{
    if (xData.size() != yData.size())
    {
        throw std::invalid_argument(
            "x and y must have the same size.");
    }

    if (xData.empty())
    {
        throw std::invalid_argument(
            "Polynomial data cannot be empty.");
    }
}

double PolynomialObjective::evaluate(
    const std::vector<double>& theta) const
{
    double loss = 0.0;

    for (size_t j = 0; j < xData.size(); ++j)
    {
        double prediction = 0.0;

        for (size_t k = 0; k < theta.size(); ++k)
        {
            prediction +=
                theta[k] * std::pow(xData[j], static_cast<double>(k));
        }

        const double error = prediction - yData[j];

        loss += error * error;
    }

    return loss / (2.0 * xData.size());
}

double PolynomialObjective::polynomial(
    const std::vector<double>& theta,
    double x) const
{
    double result = 0.0;

    for (size_t k = 0; k < theta.size(); ++k)
    {
        result += theta[k] * std::pow(x, static_cast<double>(k));
    }

    return result;
}

void PolynomialObjective::printPolynomial(const std::vector<double>& theta) const
{
    std::cout << "P(x) = ";
    bool first = true;
    for (size_t k = 0; k < theta.size(); ++k)
    {
        if (std::fabs(theta[k]) < 5e-7)
            continue;
        if (!first)
            std::cout << (theta[k] >= 0.0 ? " + " : " - ");
        else if (theta[k] < 0.0)
            std::cout << "-";
        std::cout << std::fixed << std::setprecision(6) << std::fabs(theta[k]);
        if (k >= 1)
            std::cout << "x";
        if (k >= 2)
            std::cout << "^" << k;
        first = false;
    }
    if (first)
        std::cout << "0";
    std::cout << std::endl;
}

std::vector<double> Polynomial_Approximator_Gradient_Descent(
    const std::vector<double>& x,
    const std::vector<double>& y)
{
    if (x.empty() || x.size() != y.size())
        throw std::invalid_argument("x and y must have the same non-zero size.");

    const double tolerance = 1e-5;
    GradientDescentOptions options;
    options.gradientTolerance = 1e-8;
    options.objectiveTolerance = 1e-14;
    options.stepTolerance = 1e-10;
    options.maxIterations = 10000;
    GradientDescent optimizer(options);

    double yNormSquared = 0.0;
    for (double yi : y)
        yNormSquared += yi * yi;
    const double denominator = yNormSquared == 0.0 ? 1.0 : std::sqrt(yNormSquared);

    double bestResidual = std::numeric_limits<double>::max();
    double totalCandidateTimeMs = 0.0;
    int totalIterations = 0;
    std::vector<double> bestCoefficients;
    std::vector<double> previousCoefficients;

    for (int degree = 0;
         degree < static_cast<int>(x.size()) && degree <= 50;
         ++degree)
    {
        const auto candidateStart = std::chrono::steady_clock::now();
        PolynomialObjective objective(x, y);
        std::vector<double> initialPoint(static_cast<std::size_t>(degree) + 1, 0.0);
        for (std::size_t i = 0; i < previousCoefficients.size() && i < initialPoint.size(); ++i)
            initialPoint[i] = previousCoefficients[i];
        const auto optimizerStart = std::chrono::steady_clock::now();
        OptimizationResult result = optimizer.minimize(objective, initialPoint);
        const auto optimizerEnd = std::chrono::steady_clock::now();
        if (result.theta.size() != initialPoint.size() ||
            !std::isfinite(result.objectiveValue))
            throw std::runtime_error("Gradient descent could not produce a finite polynomial fit.");

        const double residual = std::sqrt(2.0 * x.size() * result.objectiveValue) /
                                denominator;
        std::cout << "Degree = " << degree
                  << "  Residual Error = " << std::setprecision(9)
                  << residual << std::endl;
        objective.printPolynomial(result.theta);
        const auto candidateEnd = std::chrono::steady_clock::now();
        const double candidateTimeMs = std::chrono::duration<double, std::milli>(
            candidateEnd - candidateStart).count();
        totalCandidateTimeMs += std::chrono::duration<double, std::milli>(
            optimizerEnd - optimizerStart).count();
        totalIterations += result.iterations;
        std::cout << "Candidate calculation time: "
                  << candidateTimeMs
                  << " ms" << std::endl;
        std::cout << std::endl;
        previousCoefficients = result.theta;

        if (residual < bestResidual)
        {
            bestResidual = residual;
            bestCoefficients = result.theta;
        }
        if (residual <= tolerance)
            break;
    }

    const double averageTimePerIterationMs = totalIterations > 0
        ? totalCandidateTimeMs / totalIterations : 0.0;
    std::cout << "Average time per iteration: " << averageTimePerIterationMs
              << " ms (" << totalIterations << " optimizer iterations)" << std::endl;

    double maximumError = 0.0;
    PolynomialObjective bestObjective(x, y);
    for (std::size_t i = 0; i < x.size(); ++i)
        maximumError = std::max(maximumError,
            std::fabs(bestObjective.polynomial(bestCoefficients, x[i]) - y[i]));

    std::cout << std::endl
              << "Best Degree    : " << bestCoefficients.size() - 1 << std::endl
              << "Residual Error : " << std::setprecision(9) << bestResidual << std::endl
              << "Maximum Error  : " << std::setprecision(9) << maximumError << std::endl
              << "Best Fit Polynomial" << std::endl;
    bestObjective.printPolynomial(bestCoefficients);
    std::ofstream resultFile("best_polynomial_results.txt", std::ios::app);
    if (!resultFile)
        throw std::runtime_error("Unable to open best_polynomial_results.txt.");
    resultFile << std::setprecision(17)
               << "Algorithm: Gradient Descent-based\n"
               << "Best degree: " << bestCoefficients.size() - 1 << "\n"
               << "Average time per iteration (ms): " << averageTimePerIterationMs << "\n"
               << "Coefficients (ascending powers):";
    for (double coefficient : bestCoefficients)
        resultFile << " " << coefficient;
    resultFile << "\n\n";
    return bestCoefficients;
}
