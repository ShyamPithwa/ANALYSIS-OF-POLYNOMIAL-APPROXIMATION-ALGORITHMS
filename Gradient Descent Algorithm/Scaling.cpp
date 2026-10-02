#include "Scaling.h"
#include <cmath>
using namespace std;
Scaling Scaling::identity(unsigned dimension)
{
    Scaling s;
    s.scale.assign(dimension, 1.0);
    s.center.assign(dimension, 0.0);
    return s;
}

Scaling Scaling::fromBounds(const std::vector<double> &lower,
                            const std::vector<double> &upper)
{
    if (lower.size() != upper.size())
    {
        throw std::invalid_argument("Scaling bounds must have the same dimension.");
    }

    Scaling s;
    s.scale.resize(lower.size());
    s.center.resize(lower.size());

    for (size_t i = 0; i < lower.size(); ++i)
    {
        if (upper[i] <= lower[i])
        {
            throw std::invalid_argument("Upper bound must exceed lower bound for each coordinate.");
        }
        s.center[i] = 0.5 * (lower[i] + upper[i]);
        s.scale[i] = 0.5 * (upper[i] - lower[i]);
    }

    return s;
}

vector<double> Scaling::toScaled(const std::vector<double> &x) const
{
    if (x.size() != scale.size())
    {
        throw std::invalid_argument("Dimension mismatch in toScaled().");
    }

    std::vector<double> z(x.size());
    for (size_t i = 0; i < x.size(); ++i)
    {
        z[i] = (x[i] - center[i]) / scale[i];
    }
    return z;
}

std::vector<double> Scaling::toOriginal(const std::vector<double> &z) const
{
    if (z.size() != scale.size())
    {
        throw std::invalid_argument("Dimension mismatch in toOriginal().");
    }

    std::vector<double> x(z.size());
    for (size_t i = 0; i < z.size(); ++i)
    {
        x[i] = scale[i] * z[i] + center[i];
    }
    return x;
}
