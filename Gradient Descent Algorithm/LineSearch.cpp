#include "LineSearch.h"
#include "NLA_FUNCTIONS.h"
#include <cmath>

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
    int maxBacktrackingIterations)
{
    LineSearchResult result;
    double alpha = alpha0;
    const double directionalDerivative = Dot_product(g, d);

    for (int step = 0; step < maxBacktrackingIterations; ++step)
    {
        std::vector<double> candidate = Vec_Add(z, VS_Multiplication(d, alpha));

        for (double value : candidate)
        {
            if (!std::isfinite(value))
            {
                result.status = TerminationStatus::INVALID_OBJECTIVE;
                result.alpha = alpha;
                result.backtrackingSteps = step;
                return result;
            }
        }

        const double Fcandidate = F(candidate);
        if (!std::isfinite(Fcandidate))
        {
            result.status = TerminationStatus::INVALID_OBJECTIVE;
            result.alpha = alpha;
            result.backtrackingSteps = step;
            return result;
        }

        if (Fcandidate <= Fz + armijoC * alpha * directionalDerivative)
        {
            result.alpha = alpha;
            result.backtrackingSteps = step;
            result.success = true;
            result.status = TerminationStatus::CONVERGED_GRADIENT;
            return result;
        }

        alpha *= rho;
        if (alpha < minAlpha)
        {
            result.alpha = alpha;
            result.backtrackingSteps = step + 1;
            result.status = TerminationStatus::LINE_SEARCH_FAILURE;
            return result;
        }
    }

    result.alpha = alpha;
    result.backtrackingSteps = maxBacktrackingIterations;
    result.status = TerminationStatus::LINE_SEARCH_FAILURE;
    return result;
}
