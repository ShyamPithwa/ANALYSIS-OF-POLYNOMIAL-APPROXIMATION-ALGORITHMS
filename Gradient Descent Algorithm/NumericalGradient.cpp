#include "NumericalGradient.h"
#include <cmath>

NumericalGradient::NumericalGradient(Evaluator evaluator, double baseStep)
    : F(std::move(evaluator)), h(baseStep)
{
}

std::vector<double> NumericalGradient::gradient(const std::vector<double> &z) const
{
    const size_t n = z.size();
    std::vector<double> g(n, 0.0);

    for (size_t i = 0; i < n; ++i)
    {
        const double hi = h * std::max(1.0, fabs(z[i]));

        std::vector<double> zPlus = z;
        std::vector<double> zMinus = z;
        zPlus[i] += hi;
        zMinus[i] -= hi;

        const double fPlus = F(zPlus);
        const double fMinus = F(zMinus);
        g[i] = (fPlus - fMinus) / (2.0 * hi);
    }

    return g;
}
