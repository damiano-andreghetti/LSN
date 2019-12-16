#ifndef __GA__
#define __GA__
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include "random.h"

//Genetic Algorithm class
using namespace std;

class GA{
	private:
	int Np;
	int Ng;
	int Nm = 5; //5 types of mutation (implement two other)
	Random rnd;
	vector<vector<int> > Pop;
	vector<vector<double> > cities;
	vector<double> Length;
	vector<double> stats;
	vector<double> avestats;
	protected:
	
	public:
	GA(int, int);
	~GA();
	void GenerateCities(int, int);
	void GeneratePop();
	void EvaluatePop();
	void SortPop();
	int GetNp(){return Np;};
	int GetNg(){return Ng;};
	int Selection();
	int GenPBC(int);
	void PopMutation(double);
	void Mutate(int, int);
	void NewGeneration();
	bool IsInVector(int, vector<int>);
	void SaveBest(int);
	void Stats();
	void SaveStats(int);
};

#endif
