#pragma once
#include "Config.h"

class CoeffHelper
{
public:
    static double FirstDiff(const TimeIterationData& time_data, int polinomIndex);
    static double SecondDiff(const TimeIterationData& time_data, int polinomIndex);

};
