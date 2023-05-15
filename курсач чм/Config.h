#pragma once

#include <vector>
using namespace std;

enum boundary_conditions
{
	First, Second, Third
};

// struct TimeInfo
// {
// 	double StartTime = 0.0;
// 	double TimesNum = 121;
// 	double InitialStep = (double)3 / 120;
// 	double StepMultiplier = 1.2;
// };

struct TimeInfo1
{
	double StartTime = 0.0;
	double TimesNum = 61;
	double InitialStep = (double)3 / 60;
	double StepMultiplier = 1.1;
};

struct TimeIterationData
{
	TimeIterationData(double timeNow, double timeBefore1, double timeBefore2, double timeBefore3)
	{
		Times = {timeNow, timeBefore1, timeBefore2, timeBefore3};
	}

	TimeIterationData(){}

	double Delta(int diff, int startDiff = 0) const
	{
		if (diff < Times.size() && startDiff < Times.size() && diff >= 0 && startDiff >= 0)
		{
			return Times[startDiff] - Times[diff];
		}

		return -1;
	}

	vector<double> Times;

	void NextTime(double timeNow)
	{
		for (int i = Times.size() - 1; i > 0; i--)
		{
			Times[i] = Times[i - 1];
		}
		Times[0] = timeNow;
	}

	double Now() const
	{
		return Times[0];
	}
};

class Cfg
{
public:
	typedef double type;
	static const int Split = 2;
	static const int Nmat = 1;
	static const int xElem = 1 * Split;
	static const int yElem = 1 * Split;
	static constexpr type x1 = 1, x2 = 4, y1 = 1, y2 = 4;

	static constexpr type coefGEven = 2.5;
	static constexpr type coefGOdd = (double)5 / 3;

	static constexpr type coefMEven = 2.5;
	static constexpr type coefMOdd = (double)5 / 3;
};