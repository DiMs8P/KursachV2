#pragma once
class FunctionBclass
{
public:
	static double FunctionB(double r, double z, double t)
	{
		return 3*t*t + 6*t;
	}

	static double U(double r, double z, double t)
	{
		return z + t*t*t;
	}

	static double SecondBounderCondition(double r, double z, double t)
	{
		return 1;
	}
};