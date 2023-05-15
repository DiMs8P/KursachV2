#pragma once
#include <vector>
#include <iostream>
#include "Grid.h"
#include "functionB.h"

class VectorB
{
public:
	VectorB(Grid& _Grid, const vector<vector<double>>& solutions, const TimeIterationData& time_data)
	{
		this->solutions = solutions;
		this->time_data = time_data;
		GIGAGRID = &_Grid;
		CalcB();
	}

	void CalcF(int i, int k)
	{
		vector<double> first(B.size());
		for (int q = 0; q < first.size(); q++)
		{
			first[q] = solutions[solutions.size()-1][GIGAGRID->Elems[k].NodeIndex[q]] * (GIGAGRID->Mats[GIGAGRID->Elems[k].MatIndex].sigma * CoeffHelper::FirstDiff(time_data, 1) + GIGAGRID->Mats[GIGAGRID->Elems[k].MatIndex].phi * CoeffHelper::SecondDiff(time_data, 1));
		}

		vector<double> second(B.size());
		for (int q = 0; q < second.size(); q++)
		{
			second[q] = solutions[solutions.size()-2][GIGAGRID->Elems[k].NodeIndex[q]] * (GIGAGRID->Mats[GIGAGRID->Elems[k].MatIndex].sigma * CoeffHelper::FirstDiff(time_data, 2) + GIGAGRID->Mats[GIGAGRID->Elems[k].MatIndex].phi * CoeffHelper::SecondDiff(time_data, 2));
		}

		vector<double> third(B.size());
		for (int q = 0; q < third.size(); q++)
		{
			third[q] = solutions[solutions.size()-3][GIGAGRID->Elems[k].NodeIndex[q]] * (GIGAGRID->Mats[GIGAGRID->Elems[k].MatIndex].sigma * CoeffHelper::FirstDiff(time_data, 3) + GIGAGRID->Mats[GIGAGRID->Elems[k].MatIndex].phi * CoeffHelper::SecondDiff(time_data, 3));
		}
		
		for (int j = 0; j < 3; j++)
		{
			F[j] = FunctionBclass::FunctionB(GIGAGRID->Nodes[GIGAGRID->Elems[i].NodeIndex[j]].r, GIGAGRID->Nodes[GIGAGRID->Elems[i].NodeIndex[j]].z,time_data.Now()) - (first[j] + second[j] + (third[j]));
		}
	}

	void CalcB()
	{
		B.resize(3);
		
		for (int k = 0; k < GIGAGRID->Elems.size(); k++)
		{
			for (int i = 0; i < 3; i++)
			{
				B[i] = 0;
			}
			CalcF(k, k);
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3 && i != j; j++)
				{
					B[i] += F[j] * GIGAGRID->Elems[k].MatrixMMMfinale[i][j];
					B[j] += F[i] * GIGAGRID->Elems[k].MatrixMMMfinale[i][j];
				}
				B[i] += F[i] * GIGAGRID->Elems[k].MatrixMMMfinale[i][i];
			}
			
			
			GIGAGRID->Elems[k].VectorBBBfinale = B;
		}
	}
private:
	vector<vector<double>> solutions;
	TimeIterationData time_data;
	vector<double> B;
	Grid* GIGAGRID;
	double F[3];
};