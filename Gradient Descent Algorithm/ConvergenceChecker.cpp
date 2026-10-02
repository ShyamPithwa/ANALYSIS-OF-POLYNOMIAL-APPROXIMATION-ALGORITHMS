#include "ConvergenceChecker.h"
#include "NLA_FUNCTIONS.h"
#include <cmath>

TerminationStatus checkGradientConvergence(double gradientNorm,double gradientTolerance)
{
    if (gradientNorm <= gradientTolerance)
    {
        return TerminationStatus::CONVERGED_GRADIENT;
    }
    return TerminationStatus::MAX_ITERATIONS;
}

TerminationStatus checkStepConvergence(double stepNorm,double stepTolerance,const std::vector<double> &z)
{
    if (stepNorm <= stepTolerance * (1.0 + tnorm(z)))
    {
        return TerminationStatus::CONVERGED_STEP;
    }
    return TerminationStatus::MAX_ITERATIONS;
}

TerminationStatus checkObjectiveConvergence(double currentObjective,double previousObjective, double objectiveTolerance)
{
    if (std::abs(currentObjective - previousObjective) <=
        objectiveTolerance * (1.0 + std::abs(previousObjective)))
    {
        return TerminationStatus::CONVERGED_OBJECTIVE;
    }
    return TerminationStatus::MAX_ITERATIONS;
}

bool isFiniteVector(const std::vector<double> &v)
{
    for (double value : v)
    {
        if (!std::isfinite(value))
        {
            return false;
        }
    }
    return true;
}
