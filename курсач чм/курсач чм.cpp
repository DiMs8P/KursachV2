#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <set>
#include <math.h>
#include "Config.h"
#include "ElemGenerator.h"
#include "MatGenerator.h"
#include "NodeGenerator.h"
#include "GenD.h"
#include "BasicFunction.h"
#include "Grid.h"
#include "MatrixGM.h"
#include "VectorB.h"
#include "GlobalMatrixes.h"
#include "BC.h"
#include "MSG.h"
#include "Solve.h"
#include "TimeParser.h"

using namespace std;
typedef double type;
#define TIMINGSCHEME 4

int Nnode, NElem, NMat, Nbc2 = -1, Nbc3 = -1, Nbc1 = -1;
Grid _Grid;
vector<int> ig, jg;
GenD _GenD;
BC _BC;
Solve _Solve;



void Gen()
{
	GM BLIAT;
	BLIAT.GenMat();
	GE BLIAT1;
	BLIAT1.GenElem();
	GN BLIAT2;
	BLIAT2.GenNode();
}

void inputbc(BC& _BC, int& Nbc2, int& Nbc3, int& Nbc1)
{
	ifstream fin2("bc2.txt");
	ifstream fin3("bc3.txt");
	ifstream fin1("bc1.txt");
	fin2 >> Nbc2;
	_BC.BC2vector.resize(Nbc2);
	for (int i = 0; i < Nbc2; i++)
	{
		fin2 >> _BC.BC2vector[i].versh[0] >> _BC.BC2vector[i].versh[1] >> _BC.BC2vector[i].theta[0] >> _BC.BC2vector[i].theta[1];
	}
	fin3 >> Nbc3;
	_BC.BC3vector.resize(Nbc3);
	for (int i = 0; i < Nbc3; i++)
	{
		fin3 >> _BC.BC3vector[i].versh[0] >> _BC.BC3vector[i].versh[1] >> _BC.BC3vector[i].Beta >> _BC.BC3vector[i].U[0] >> _BC.BC3vector[i].U[1];
	}
	fin1 >> Nbc1;
	_BC.BC1vector.resize(Nbc1);
	for (int i = 0; i < Nbc1; i++)
	{
		fin1 >> _BC.BC1vector[i].versh[0] >> _BC.BC1vector[i].versh[1] >> _BC.BC1vector[i].U;
	}
}

void InputGrid(Grid& _Grid, int& Nnode, int& NElem, int& NMat)
{
	ifstream finN("Node.txt");
	ifstream finE("Elem.txt");
	ifstream finM("Mat.txt");
	finN >> Nnode;
	_Grid.Nodes.resize(Nnode);
	for (size_t i = 0; i < Nnode; i++)
	{
		finN >> _Grid.Nodes[i].r >> _Grid.Nodes[i].z;
	}
	finE >> NElem;
	_Grid.Elems.resize(NElem);
	for (size_t i = 0; i < NElem; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			finE >> _Grid.Elems[i].NodeIndex[j];
		}
		finE >> _Grid.Elems[i].MatIndex;
	}
	finM >> NMat;
	_Grid.Mats.resize(NMat);
	for (size_t i = 0; i < NMat; i++)
	{
		finM >> _Grid.Mats[i].gamma >> _Grid.Mats[i].sigma >> _Grid.Mats[i].phi >> _Grid.Mats[i].L;
	}
	_Grid.height = _Grid.Nodes[_Grid.Elems[0].NodeIndex[2]].z - _Grid.Nodes[0].z;
	_Grid.width = _Grid.Nodes[1].r - _Grid.Nodes[0].r;
}

void InitializeSolutions(const vector<double>& times, vector<vector<double>>& solutions)
{
	for(int i = 0; i < TIMINGSCHEME - 1; i++)
	{
		for (int j = 0; j < solutions[i].size(); j++)
		{
			solutions[i][j] = FunctionBclass::U(_Grid.Nodes[j].r, _Grid.Nodes[j].z, times[i]);
		}
	}
}

void UpdateFirstBoundaryCondition(const TimeIterationData& time_info)
{
	for (int i = 0; i < _BC.BC1vector.size(); i++)
	{
		_BC.BC1vector[i].U = FunctionBclass::U(_Grid.Nodes[_BC.BC1vector[i].versh[0]].r, _Grid.Nodes[_BC.BC1vector[i].versh[0]].z, time_info.Now());
	}
}

void UpdateSecondBoundaryCondition(const TimeIterationData& time_info)
{
	// for (int i = 0; i < _BC.BC1vector.size(); i++)
	// {
	// 	_BC.BC2vector[i].theta[0] = FunctionBclass::SecondBounderCondition(_Grid.Nodes[_BC.BC1vector[i].versh[0]].r, _Grid.Nodes[_BC.BC1vector[i].versh[0]].z, time_info.Now());
	// 	_BC.BC2vector[i].theta[1] = FunctionBclass::SecondBounderCondition(_Grid.Nodes[_BC.BC1vector[i].versh[0]].r, _Grid.Nodes[_BC.BC1vector[i].versh[0]].z, time_info.Now());
	// }
}

void UpdateThirdBoundaryCondition(const TimeIterationData& time_info)
{
	
}

void UpdateBoundaryConditions(const TimeIterationData& time_info, const vector<boundary_conditions>& boundary_conditions)
{
	for (int i = 0; i < boundary_conditions.size(); i++)
	{
	    switch (boundary_conditions[i])
	    {
	    case First:
	        {
	            UpdateFirstBoundaryCondition(time_info);
	            break;
	        }
	
	    case Second:
	        {
	            UpdateSecondBoundaryCondition(time_info);
	            break;
	        }
	
	    case Third:
	        {
	            UpdateThirdBoundaryCondition(time_info);
	            break;
	        }
	    }
	}
}

void PrintAt(double r, double z, const vector<double>& solutions)
{
	double soooooolve;
	double koefdiag, koefdot;
	int koefZ, koefR, IndexElem;
	koefZ = z / _Grid.height;
	koefR = r / _Grid.width;
	if (koefR > (Cfg::xElem - 1))
		koefR = koefR - 1;
	if (koefZ > (Cfg::yElem - 1))
		koefZ = koefZ - 1;
	koefdot = ((koefR + 1) * _Grid.width - (r - Cfg::x1)) / ((z - Cfg::y1) - koefZ * _Grid.height);
	if (z - Cfg::y1 - koefZ * _Grid.height <= 0.0001)
		koefdot = 1;

	koefdiag = _Grid.width / _Grid.height;

	if (koefdot > koefdiag)
		IndexElem = ((koefZ * Cfg::xElem) + koefR) * 2;
	else
		IndexElem = (((koefZ * Cfg::xElem) + koefR) * 2) + 1;
	soooooolve = 0;
	for (int i = 0; i < 3; i++)
	{
		soooooolve += solutions[_Grid.Elems[IndexElem].NodeIndex[i]] * _Grid.Elems[IndexElem].BF[i].functionIn(r, z);
	}
	cout << setprecision(15) << soooooolve << std::endl;
}

int main()
{
	cout.imbue(std::locale("German_germany"));
	 Gen();
	 InputGrid(_Grid, Nnode, NElem, NMat);
	 _GenD.D(_Grid, NElem);
	
	vector<vector<double>> solutions(TIMINGSCHEME - 1, std::vector<double>(_Grid.Nodes.size()));
	
	TimeInfo1 time_info;
	TimeParser time_parser(time_info);
	vector<double> times;
	time_parser.Parse(times);
	
	InitializeSolutions(times, solutions);
	
	TimeIterationData time_iteration_data(times[2], times[1], times[0], 0);
	MatrixMG _MatrixMG(_Grid);
	
	inputbc(_BC, Nbc2, Nbc3, Nbc1);
	cout << setprecision(15) << "U" << std::endl;
	for (int i = TIMINGSCHEME - 1; i < time_info.TimesNum; i++)
	{
		time_iteration_data.NextTime(times[i]);
		UpdateBoundaryConditions(time_iteration_data, {First, Second, Third});
		
		VectorB _VectorB(_Grid, solutions, time_iteration_data);
		GlobalMatrix _GlobalMatrix(_Grid, NElem, Nnode, time_iteration_data);
		
		_BC.primeniaemKraevble(_Grid, _GlobalMatrix, Nbc2, Nbc3, Nbc1, Nnode);
		SLAU _SLAU(_GlobalMatrix, Nnode);
		solutions.push_back(_SLAU.q);
		//_Solve.finale(_Grid, _SLAU);

		if (i == 20 || i == 40 || i == 60)
		{
			cout << "t = " << times[i] << endl;
			PrintAt(2.5, 2.5, solutions[i]);
			PrintAt(2.5, 4, solutions[i]);
		}
	}

	cout << endl;
	cout << setprecision(15) << "U*" << std::endl;
	
	for (int i = 1; i <= 3; i++)
	{
		cout << "t = " << times[i * 20] << endl;
		cout << setprecision(15) << FunctionBclass::U(2.5, 2.5, times[i * 20]) << std::endl;
		cout << setprecision(15) << FunctionBclass::U(2.5, 4, times[i * 20]) << std::endl;
	}

	
	return 0;
}


