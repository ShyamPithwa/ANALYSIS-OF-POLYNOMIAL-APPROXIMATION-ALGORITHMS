#pragma once
#ifndef GRADIENT_DESCENT_OPTIONS_H
#define GRADIENT_DESCENT_OPTIONS_H

#include <vector>

struct GradientDescentOptions
{
    bool useScaling = false;
    std::vector<double> lowerBounds;
    std::vector<double> upperBounds;

    double baseFiniteDifferenceStep = 1e-5;
    double initialAlpha = 1.0;
    bool warmStartAlpha = true;
    double backtrackingRho = 0.5;
    double armijoC = 1e-4;

    double gradientTolerance = 1e-8;
    double objectiveTolerance = 1e-12;
    double stepTolerance = 1e-10;
    int maxIterations = 10000;

    double minAlpha = 1e-16;
    int maxBacktrackingIterations = 100;
    bool storeIterationHistory = false;
};

#endif
