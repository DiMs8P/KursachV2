#include "CoeffHelper.h"

double CoeffHelper::FirstDiff(const TimeIterationData& time_data, int polinomIndex)
{
    switch (polinomIndex)
    {
    case 0:
        {
            return (time_data.Delta(2) * time_data.Delta(1) + time_data.Delta(3) * time_data.Delta(1) + time_data.Delta(2) * time_data.Delta(3))
        / (time_data.Delta(1) * time_data.Delta(2) * time_data.Delta(3));
        }

    case 1:
        {
            return (time_data.Delta(2) * time_data.Delta(3)) / (time_data.Delta(3, 1) * time_data.Delta(2, 1) * time_data.Delta(0, 1));
        }

    case 2:
        {
            return (time_data.Delta(1) * time_data.Delta(3))/ (time_data.Delta(3, 2) * time_data.Delta(1, 2) * time_data.Delta(0, 2));
        }

    case 3:
        {
            return (time_data.Delta(1) * time_data.Delta(2))/ (time_data.Delta(2, 3) * time_data.Delta(1, 3) * time_data.Delta(0, 3));
        }
    }
}

double CoeffHelper::SecondDiff(const TimeIterationData& time_data, int polinomIndex)
{
    switch (polinomIndex)
    {
    case 0:
        {
            return (2 * (time_data.Delta(1) + time_data.Delta(2) + time_data.Delta(3)))
            / (time_data.Delta(1) * time_data.Delta(2) * time_data.Delta(3));
        }

    case 1:
        {
            return (2 * (time_data.Delta(2) + time_data.Delta(3)))
            / (time_data.Delta(3, 1) * time_data.Delta(2, 1) * time_data.Delta(0, 1));
        }

    case 2:
        {
            return (2 * (time_data.Delta(1) + time_data.Delta(3)))
            / (time_data.Delta(3, 2) * time_data.Delta(1, 2) * time_data.Delta(0, 2));
        }

    case 3:
        {
            return (2 * (time_data.Delta(2) + time_data.Delta(1)))
            / (time_data.Delta(2, 3) * time_data.Delta(1, 3) * time_data.Delta(0, 3));
        }
    }
}
