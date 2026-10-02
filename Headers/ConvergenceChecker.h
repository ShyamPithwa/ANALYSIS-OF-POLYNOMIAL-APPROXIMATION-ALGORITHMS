#pragma once
#ifndef CONVERGENCE_CHECKER_H
#define CONVERGENCE_CHECKER_H

#include "OptimizationResult.h"
#include <vector>

struct ConvergenceThresholds
{
    double gradientTolerance = 1e-8;
    double objectiveTolerance = 1e-12;
    double stepTolerance = 1e-10;
};

TerminationStatus checkGradientConvergence(double gradientNorm,
                                           double gradientTolerance);

TerminationStatus checkStepConvergence(double stepNorm,
                                       double stepTolerance,
                                       const std::vector<double> &z);

TerminationStatus checkObjectiveConvergence(double currentObjective,
                                            double previousObjective,
                                            double objectiveTolerance);

bool isFiniteVector(const std::vector<double> &v);

#endif
