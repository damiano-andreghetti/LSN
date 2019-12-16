#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "blocchi.h"
#include "GA.h"
#include <vector>
#define _USE_MATH_DEFINES 
using namespace std;

/*
CREARE CLASSE PER ALGORITMO GENETICO
-cosi uso un unico oggetto random
*/	

/*
int Nc = 30; //# of cities
int Np = 900; //# of members in the poulation
vector<vector<double>> cities(Nc);
vector<vector<int>> Pop(Np);
vector<double> Length(Np);

int Selection(int Np){
	Random rnd;
	int seed[4];
  int p1, p2;
  ifstream Primes("Primes");
  if (Primes.is_open()){
     Primes >> p1 >> p2 ;
  } else cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();
  ifstream input("seed.in");
  string property;
  if (input.is_open()){
     while ( !input.eof() ){
        input >> property;
        if( property == "RANDOMSEED" ){
           input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
           rnd.SetRandom(seed,p1,p2);
        }
     }
     input.close();
  } else cerr << "PROBLEM: Unable to open seed.in" << endl;
	es = 0.5;
	int j = Np* int(pow(rnd.Rannyu(), es));
	rnd.SaveSeed();
	return j;
}


void GenerateCities(int N){
	Random rnd;
	int seed[4];
  int p1, p2;
  ifstream Primes("Primes");
  if (Primes.is_open()){
     Primes >> p1 >> p2 ;
  } else cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();
  ifstream input("seed.in");
  string property;
  if (input.is_open()){
     while ( !input.eof() ){
        input >> property;
        if( property == "RANDOMSEED" ){
           input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
           rnd.SetRandom(seed,p1,p2);
        }
     }
     input.close();
  } else cerr << "PROBLEM: Unable to open seed.in" << endl;
	for(int i = 0; i < N; i++){
		double teta = rnd.Rannyu(0., 2.*M_PI);
		double xc = cos(teta);
		double yc = sin(teta);
		//cout << xc << "," << yc << endl;
		cities[i] = vector<double>(2);
		cities[i][0]=xc;
		cities[i][1]=yc;
	}
	rnd.SaveSeed();
}

void GeneratePopulation(int Np, int Ng){
	Random rnd;
	int seed[4];
  int p1, p2;
  ifstream Primes("Primes");
  if (Primes.is_open()){
     Primes >> p1 >> p2 ;
  } else cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();
  ifstream input("seed.in");
  string property;
  if (input.is_open()){
     while ( !input.eof() ){
        input >> property;
        if( property == "RANDOMSEED" ){
           input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
           rnd.SetRandom(seed,p1,p2);
        }
     }
     input.close();
  } else cerr << "PROBLEM: Unable to open seed.in" << endl;
	for(int i=0; i <Np; i++){
		Pop[i] = vector<int> (Ng);
		Pop[i][0]=int(Ng*rnd.Rannyu());
		for(int j = 1; j < Ng; j++){
			int a;
			bool already_used;
			do{
				already_used=0;
				a = int(Ng*rnd.Rannyu());
				for(int k = 0; k <j; k++) if(Pop[i][k]==a) already_used=1;
			}
			while(already_used);
			Pop[i][j] = a;
		}
		//cout << "generated individual n: " << i+1 << endl;
	}
	rnd.SaveSeed();
}

void EvaluatePop(int Np, int Ng){
	//prepare length vector
	for(int i = 0; i < Np; i++){
			double l = 0;
			for(int j = 0; j < Ng; j++){
				if(j==Ng-1) l+= distance(cities[Pop[i][j]][0],cities[Pop[i][j]][1],cities[Pop[i][0]][0],cities[Pop[i][0]][1]);
				else l+= distance(cities[Pop[i][j]][0],cities[Pop[i][j]][1],cities[Pop[i][j+1]][0],cities[Pop[i][j+1]][1]);
			}
			Length[i]=l;	
			//cout << l << endl;
		} 
}

void SortPop(int Np){
	bool sorted;
	do{
		sorted=1;
		for(int i =1; i < Np; i++){
			if(Length[i]<Length[i-1]){ //swap elements
				sorted=0;
				double v = Length[i];
				Length[i]=Length[i-1];
				Length[i-1]=v;
				Pop[i].swap(Pop[i-1]);
			}
		}
	}
	while(!sorted);
}
*/

int main(){
	GA ga(900,30);
	cout << "inizializzato GA" << endl;
	ga.GenerateCities(30,0);
	cout << "cities on unitary circle generated" << endl;
	ga.GeneratePop();
	cout << "population generated" << endl;
	ga.EvaluatePop();
	cout << "population evaluated" << endl;
	ga.SortPop();
	cout << "populaton sorted" << endl;		
	ga.Stats();
	for(int i = 0; i < 100; i++){
		cout << "creating generation no." << i+1 << endl;	
		ga.NewGeneration();
		ga.PopMutation(0.1);
		ga.EvaluatePop();
		ga.SortPop();
		ga.Stats();
	}
	ga.SaveBest(1);
	ga.SaveStats(1);
	GA ga2(900,30);
	cout << "inizializzato GA" << endl;
	ga2.GenerateCities(30,1);
	cout << "cities in a square generated" << endl;
	ga2.GeneratePop();
	cout << "population generated" << endl;
	ga2.EvaluatePop();
	cout << "population evaluated" << endl;
	ga2.SortPop();
	cout << "populaton sorted" << endl;		
	ga2.Stats();
	for(int i = 0; i < 100; i++){
		cout << "creating generation no." << i+1 << endl;	
		ga2.NewGeneration();
		ga2.PopMutation(0.1);
		ga2.EvaluatePop();
		ga2.SortPop();
		ga2.Stats();
	}
	ga2.SaveBest(2);
	ga2.SaveStats(2);
	return 0;
}


/*
unsigned long int N = 200000;
unsigned long int M = 20;
unsigned long int L = N/M;
double xn;
double x0=0;
double delta = 2.95; //delta per gli step per metropolis
double xn1;
double a0 = 1.;

double potential(double x){
	double v = pow(x, 4.) - 2.5*x*x;
	return v;
}

double MODPSI(double x, double y, double z){ //funzione d'onda atomo idrogeno, nlm=100
	double u = pow(a0, -3.)*exp(-2.*sqrt(x*x+y*y+z*z)/a0)/M_PI;
	return u;
}

double MODPSI1D(double m,double s, double x){
	double u = exp(-((x-m)*(x-m))/(s*s))+exp(-(x*x+m*m)/(s*s))+exp(-((x+m)*(x+m))/(s*s));
	return u;
}
double MODPSI2(double x, double y, double z){ //funzione d'onda atomo idrogeno, nlm=210
	double u = z*z*exp(-sqrt(x*x+y*y+z*z)/a0)/(M_PI*32.*a0*a0*a0*a0*a0);
	return u;
}

int main(){
	Random rnd;
	Blocchi  bl(N, M, "output.dat");
  int seed[4];
  int p1, p2;
  ifstream Primes("Primes");
  if (Primes.is_open()){
     Primes >> p1 >> p2 ;
  } else cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();
  ifstream input("seed.in");
  string property;
  if (input.is_open()){
     while ( !input.eof() ){
        input >> property;
        if( property == "RANDOMSEED" ){
           input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
           rnd.SetRandom(seed,p1,p2);
        }
     }
     input.close();
  } else cerr << "PROBLEM: Unable to open seed.in" << endl;
	
	cout << "probablità di transizione uniforme" << endl;
	int acc=0;
	//metropolis step START	
	//prima faccio 1000000 di punti a vuoto;
	ofstream pos;
	pos.open("pos.xyz");
	unsigned long int acc=0;
	for(int p=0; p < M; p++){
		//cout << "reset position, start new block" <<endl;
		xn[0]=x0[0];xn[1]=x0[1];xn[2]=x0[2];
		double a[3] = {0.,0.,0.};
		for(int q = 0; q<L; q++){
			xn1[0]=xn[0]; xn1[1]=xn[1]; xn1[2]=xn[2];
			for(int i = 0; i < 3; i++){
				xn1[i]=xn[i]+rnd.Rannyu(-delta, delta);
				a[i]=min(1.,MODPSI(xn1[0], xn1[1], xn1[2])/MODPSI(xn[0], xn[1], xn[2]));
				if(rnd.Rannyu() <= a[i]){ //accept
					acc++;
					xn[i]=xn1[i];
				}
				else { //reject
					xn1[i]=xn[i];
				} 
			}
			if(((p*L)+q)%100 )pos << (p*L)+q << " " << xn[0] << " " << xn[1] << " " << xn[2] << endl; //plot some position (1 over 10)
			bl.addThrow(sqrt(xn[0]*xn[0]+xn[1]*xn[1]+xn[2]*xn[2]));
		}
		
	}
	pos.close();
	//metropolis step END

	double mu =0.8;
	double sigma = 0.8;
	for(int p=0; p < M; p++){
		xn=x0;
		for(int q = 0; q<L; q++){
			xn1=xn+rnd.Rannyu(-delta,delta);
			double a =min(1.,MODPSI1D(mu, sigma, xn1)/MODPSI1D(mu, sigma, xn));
			if(rnd.Rannyu() <= a){
				acc++;
				xn=xn1;
			}
			else{} //reject
			double stima_H = potential(xn);
			bl.addThrow(stima_H);
			}
	}
	cout << acc*100./N << endl;
  rnd.SaveSeed();
  return 0;
}
*/

