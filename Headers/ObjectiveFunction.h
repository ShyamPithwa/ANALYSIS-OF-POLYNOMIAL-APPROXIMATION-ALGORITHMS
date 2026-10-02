#pragma once
#ifndef OBJECTIVE_FUNCTION_H
#define OBJECTIVE_FUNCTION_H

#include <vector>

class ObjectiveFunction
{
public:
    virtual ~ObjectiveFunction() = default;
    virtual double evaluate(const std::vector<double> &theta) const = 0;
};

#endif
