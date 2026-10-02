#pragma once
#ifndef OPTIMIZATION_RESULT_H
#define OPTIMIZATION_RESULT_H

#include <vector>
#include <string>

enum class TerminationStatus
{
    CONVERGED_GRADIENT,
    CONVERGED_OBJECTIVE,
    CONVERGED_STEP,
    MAX_ITERATIONS,
    LINE_SEARCH_FAILURE,
    INVALID_OBJECTIVE,
    INVALID_GRADIENT
};

inline std::string terminationStatusToString(TerminationStatus status)
{
    switch (status)
    {
    case TerminationStatus::CONVERGED_GRADIENT:
        return "CONVERGED_GRADIENT";
    case TerminationStatus::CONVERGED_OBJECTIVE:
        return "CONVERGED_OBJECTIVE";
    case TerminationStatus::CONVERGED_STEP:
        return "CONVERGED_STEP";
    case TerminationStatus::MAX_ITERATIONS:
        return "MAX_ITERATIONS";
    case TerminationStatus::LINE_SEARCH_FAILURE:
        return "LINE_SEARCH_FAILURE";
    case TerminationStatus::INVALID_OBJECTIVE:
        return "INVALID_OBJECTIVE";
    case TerminationStatus::INVALID_GRADIENT:
        return "INVALID_GRADIENT";
    default:
        return "UNKNOWN";
    }
}

struct IterationRecord
{
    int iteration;
    double objective;
    double gradientNorm;
    double alpha;
    double stepNorm;
};

struct OptimizationResult
{
    std::vector<double> theta;
    double objectiveValue = 0.0;
    double gradientNorm = 0.0;
    double stepNorm = 0.0;
    int iterations = 0;
    TerminationStatus status = TerminationStatus::MAX_ITERATIONS;
    double acceptedAlpha = 0.0;
    int backtrackingSteps = 0;
    std::vector<IterationRecord> history;
};

#endif
