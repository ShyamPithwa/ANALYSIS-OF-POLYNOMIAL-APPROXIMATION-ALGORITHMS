#pragma once
#ifndef GRADIENT_DESCENT_H
#define GRADIENT_DESCENT_H

#include "GradientDescentOptions.h"
#include "ObjectiveFunction.h"
#include "OptimizationResult.h"
#include "Scaling.h"

class GradientDescent
{
public:
    explicit GradientDescent(const GradientDescentOptions &options = GradientDescentOptions());

    OptimizationResult minimize(const ObjectiveFunction &objective,
                                const std::vector<double> &initialPoint) const;

private:
    GradientDescentOptions options;

    void finalizeResult(OptimizationResult &result,
                        const std::vector<double> &z,
                        const Scaling &scaling,
                        double objectiveValue,
                        double gradientNorm,
                        double stepNorm,
                        int iterations,
                        TerminationStatus status,
                        double acceptedAlpha,
                        int backtrackingSteps) const;
};

#endif
