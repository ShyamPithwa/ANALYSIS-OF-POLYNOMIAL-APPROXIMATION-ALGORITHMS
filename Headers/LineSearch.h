#pragma once
#ifndef LINE_SEARCH_H
#define LINE_SEARCH_H

#include "OptimizationResult.h"
#include <functional>
#include <vector>

struct LineSearchResult
{
    double alpha = 0.0;
    int backtrackingSteps = 0;
    TerminationStatus status = TerminationStatus::LINE_SEARCH_FAILURE;
    bool success = false;
};

LineSearchResult backtrackingArmijo(
    const std::function<double(const std::vector<double> &)> &F,
    const std::vector<double> &z,
    const std::vector<double> &g,
    const std::vector<double> &d,
    double Fz,
    double alpha0,
    double armijoC,
    double rho,
    double minAlpha,
    int maxBacktrackingIterations);

#endif
