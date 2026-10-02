#pragma once
#ifndef OPTIMIZATION_SCALING_H
#define OPTIMIZATION_SCALING_H

#include <vector>
#include <stdexcept>

class Scaling
{
public:
    static Scaling identity(unsigned dimension);
    static Scaling fromBounds(const std::vector<double> &lower,
                              const std::vector<double> &upper);

    std::vector<double> toScaled(const std::vector<double> &x) const;
    std::vector<double> toOriginal(const std::vector<double> &z) const;

    unsigned dimension() const { return scale.size(); }
    const std::vector<double> &getScale() const { return scale; }
    const std::vector<double> &getCenter() const { return center; }

private:
    std::vector<double> scale;
    std::vector<double> center;
};

#endif
