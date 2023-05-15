#pragma once
#include"MatrixGM.h"

class GlobalMatrix
{
public:

	void SumMatrixMG(const TimeIterationData& time_data)
		{
			double sum;

			double firstCoef = CoeffHelper::FirstDiff(time_data, 0);
			double secondCoef = CoeffHelper::SecondDiff(time_data, 0);
			for (int k = 0; k < GIGAGRID->Elems.size(); k++) 
			{
				GIGAGRID->Elems[k].GIGAMATRIX.resize(3);
				for (int i = 0; i < 3; i++)
				{
					GIGAGRID->Elems[k].GIGAMATRIX[i].resize(3);
				}
				for (int i = 0; i < 3; i++) {
					for (int j = 0; j <= i; j++) {
						GIGAGRID->Elems[k].GIGAMATRIX[i][j] = GIGAGRID->Elems[k].MatrixGGGfinale[i][j] + (GIGAGRID->Elems[k].MatrixMMMfinale[i][j]
						* GIGAGRID->Mats[GIGAGRID->Elems[k].MatIndex].sigma * firstCoef
						+ GIGAGRID->Elems[k].MatrixMMMfinale[i][j]
						* GIGAGRID->Mats[GIGAGRID->Elems[k].MatIndex].phi * secondCoef);
					}
				}
			}
	}

	GlobalMatrix(Grid& _Grid, int NElem, int Nnode, const TimeIterationData& time_data)
	{
		GIGAGRID = &_Grid;
		SumMatrixMG(time_data);
		portrait(_Grid, NElem, Nnode);
		GlobalGMtriangle.resize(ig[Nnode - 1]);
		GlobalGMdiag.resize(Nnode);
		Globalvector.resize(Nnode);
		vstavka(_Grid, NElem);
	}
	vector<double> GlobalGMtriangle, GlobalGMdiag, Globalvector;
	vector<int> ig, jg;
	Grid* GIGAGRID;
	void portrait(Grid& _Grid, int NElem, int Nnode)
	{

		versh.resize(Nnode);
		ig.resize(Nnode);

		for (int i = 0; i < NElem; i++)
		{
			versh[_Grid.Elems[i].NodeIndex[1]].insert(_Grid.Elems[i].NodeIndex[0]);
			for (int j = 0; j < 2; j++)
			{
				versh[_Grid.Elems[i].NodeIndex[2]].insert(_Grid.Elems[i].NodeIndex[j]);
			}
		}
		ig[0] = 0;
		for (int i = 1; i < Nnode; i++)
		{
			ig[i] = ig[i - 1] + versh[i].size();

		}
		for (int i = 0; i < versh.size(); i++)
		{
			for (set<int>::iterator it = versh[i].begin(); it != versh[i].end(); ++it)
				jg.push_back(*it);
		}
	}

	void vstavka(Grid _Grid, int NElem)
	{
		for (int i = 0; i < NElem; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				GlobalGMdiag[_Grid.Elems[i].NodeIndex[j]] += _Grid.Elems[i].GIGAMATRIX[j][j];
			}
		}
		for (int i = 0; i < NElem; i++)
		{
			for (int k = ig[_Grid.Elems[i].NodeIndex[1]-1]; k < ig[_Grid.Elems[i].NodeIndex[1]]; k++)
			{
				if(jg[k] == _Grid.Elems[i].NodeIndex[0])
					GlobalGMtriangle[k] += _Grid.Elems[i].GIGAMATRIX[1][0];
			}		
			for (int j = 0; j < 2; j++)
			{
				for (int k = ig[_Grid.Elems[i].NodeIndex[2] - 1]; k < ig[_Grid.Elems[i].NodeIndex[2]]; k++)
				{
					if (jg[k] == _Grid.Elems[i].NodeIndex[j]) {
						GlobalGMtriangle[k] += _Grid.Elems[i].GIGAMATRIX[2][j];
						k += 3;
					}
				}
			}
		}
		for (int i = 0; i < NElem; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				Globalvector[_Grid.Elems[i].NodeIndex[j]] += _Grid.Elems[i].VectorBBBfinale[j]; 
			}
		}
	}

private:
	vector<set<int>> versh;
};
