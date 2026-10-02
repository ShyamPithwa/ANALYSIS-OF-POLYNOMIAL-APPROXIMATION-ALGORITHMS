#include "GradientDescent.h"
#include "ConvergenceChecker.h"
#include "LineSearch.h"
#include "NumericalGradient.h"
#include "NLA_FUNCTIONS.h"
#include <cmath>

GradientDescent::GradientDescent(const GradientDescentOptions &options)
    : options(options)
{
}

void GradientDescent::finalizeResult(OptimizationResult &result,
                                     const std::vector<double> &z,
                                     const Scaling &scaling,
                                     double objectiveValue,
                                     double gradientNorm,
                                     double stepNorm,
                                     int iterations,
                                     TerminationStatus status,
                                     double acceptedAlpha,
                                     int backtrackingSteps) const
{
    result.theta = scaling.toOriginal(z);
    result.objectiveValue = objectiveValue;
    result.gradientNorm = gradientNorm;
    result.stepNorm = stepNorm;
    result.iterations = iterations;
    result.status = status;
    result.acceptedAlpha = acceptedAlpha;
    result.backtrackingSteps = backtrackingSteps;
}

OptimizationResult GradientDescent::minimize(const ObjectiveFunction &objective,
                                             const std::vector<double> &initialPoint) const
{
    OptimizationResult result;

    if (initialPoint.empty())
    {
        result.status = TerminationStatus::INVALID_OBJECTIVE;
        return result;
    }

    Scaling scaling = options.useScaling
                          ? Scaling::fromBounds(options.lowerBounds, options.upperBounds)
                          : Scaling::identity(static_cast<unsigned>(initialPoint.size()));

    if (options.useScaling &&
        (options.lowerBounds.size() != initialPoint.size() ||
         options.upperBounds.size() != initialPoint.size()))
    {
        result.status = TerminationStatus::INVALID_OBJECTIVE;
        return result;
    }

    std::vector<double> z = scaling.toScaled(initialPoint);

    auto evaluateScaled = [&](const std::vector<double> &zval) -> double
    {
        std::vector<double> x = scaling.toOriginal(zval);
        return objective.evaluate(x);
    };

    NumericalGradient gradientProvider(evaluateScaled, options.baseFiniteDifferenceStep);

    double trialAlpha = options.initialAlpha;
    double currentObjective = evaluateScaled(z);

    if (!std::isfinite(currentObjective))
    {
        result.status = TerminationStatus::INVALID_OBJECTIVE;
        result.theta = initialPoint;
        result.objectiveValue = currentObjective;
        return result;
    }

    for (int iteration = 0; iteration < options.maxIterations; ++iteration)
    {
        const double Fk = currentObjective;

        std::vector<double> gk = gradientProvider.gradient(z);
        const double gradientNorm = tnorm(gk);

        if (!isFiniteVector(gk))
        {
            finalizeResult(result, z, scaling, Fk, gradientNorm, 0.0, iteration,
                           TerminationStatus::INVALID_GRADIENT, 0.0, 0);
            return result;
        }

        if (checkGradientConvergence(gradientNorm, options.gradientTolerance) ==
            TerminationStatus::CONVERGED_GRADIENT)
        {
            finalizeResult(result, z, scaling, Fk, gradientNorm, 0.0, iteration,
                           TerminationStatus::CONVERGED_GRADIENT, 0.0, 0);
            return result;
        }

        std::vector<double> dk = VS_Multiplication(gk, -1.0);

        LineSearchResult lineSearch = backtrackingArmijo(
            evaluateScaled,
            z,
            gk,
            dk,
            Fk,
            trialAlpha,
            options.armijoC,
            options.backtrackingRho,
            options.minAlpha,
            options.maxBacktrackingIterations);

        if (!lineSearch.success)
        {
            finalizeResult(result, z, scaling, Fk, gradientNorm, 0.0, iteration,
                           lineSearch.status, lineSearch.alpha, lineSearch.backtrackingSteps);
            return result;
        }

        const double acceptedAlpha = lineSearch.alpha;
        std::vector<double> zNext = Vec_Add(z, VS_Multiplication(dk, acceptedAlpha));
        const double stepNorm = tnorm(Vec_Sub(zNext, z));

        if (!isFiniteVector(zNext))
        {
            finalizeResult(result, z, scaling, Fk, gradientNorm, stepNorm, iteration,
                           TerminationStatus::INVALID_OBJECTIVE, acceptedAlpha,
                           lineSearch.backtrackingSteps);
            return result;
        }

        const double nextObjective = evaluateScaled(zNext);
        if (!std::isfinite(nextObjective))
        {
            finalizeResult(result, z, scaling, Fk, gradientNorm, stepNorm, iteration,
                           TerminationStatus::INVALID_OBJECTIVE, acceptedAlpha,
                           lineSearch.backtrackingSteps);
            return result;
        }

        if (options.storeIterationHistory)
        {
            IterationRecord record;
            record.iteration = iteration;
            record.objective = Fk;
            record.gradientNorm = gradientNorm;
            record.alpha = acceptedAlpha;
            record.stepNorm = stepNorm;
            result.history.push_back(record);
        }

        if (checkStepConvergence(stepNorm, options.stepTolerance, z) ==
            TerminationStatus::CONVERGED_STEP)
        {
            finalizeResult(result, zNext, scaling, nextObjective, gradientNorm, stepNorm,
                           iteration + 1, TerminationStatus::CONVERGED_STEP, acceptedAlpha,
                           lineSearch.backtrackingSteps);
            return result;
        }

        if (checkObjectiveConvergence(nextObjective, Fk, options.objectiveTolerance) ==
            TerminationStatus::CONVERGED_OBJECTIVE)
        {
            finalizeResult(result, zNext, scaling, nextObjective, gradientNorm, stepNorm,
                           iteration + 1, TerminationStatus::CONVERGED_OBJECTIVE, acceptedAlpha,
                           lineSearch.backtrackingSteps);
            return result;
        }

        z = zNext;
        currentObjective = nextObjective;
        trialAlpha = options.warmStartAlpha ? acceptedAlpha : options.initialAlpha;
    }

    const std::vector<double> finalGradient = gradientProvider.gradient(z);
    finalizeResult(result, z, scaling, currentObjective, tnorm(finalGradient), 0.0,
                   options.maxIterations, TerminationStatus::MAX_ITERATIONS, 0.0, 0);
    return result;
}
