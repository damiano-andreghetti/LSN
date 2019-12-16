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
#include <utility>
#include <tuple>

//Genetic Algorithm class
using namespace std;

class SA{
	private:
	double sT, cT; //starting temperature, current temperature
	int is, icT;  //temperature step index, current temperature step index
	int NTS, NSS; //number of temperarure step, number of sub steps for each temperature
	int Ng;
	double dT;
	Random rnd;
	vector<vector<double> > cities;
	vector<int> config;
	vector<int> config1;
	vector<int> bst_config;
	double L_bst=2000.;
	//vector<double> Length; 
	//vector<double> stats;
	//vector<double> avestats;
	protected:
	
	public:
	SA(double, int, int, int, int);
	~SA();
	void GenerateCities(int);
	double Run();
	pair<double,double> evaluate();
	int GenPBC(int);
	void Mutate(int);
	void SaveBest(int);
	/*
	void GeneratePop();
	void EvaluatePop();
	void SortPop();
	int GetNp(){return Np;};
	int GetNg(){return Ng;};
	int Selection();
	void NewGeneration();
	bool IsInVector(int, vector<int>);
	void Stats();
	void SaveStats(int);
	*/
};

#endif
