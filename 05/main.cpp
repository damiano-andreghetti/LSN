#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "blocchi.h"

#define _USE_MATH_DEFINES 
using namespace std;


int N = 1000000;
int M = 100;	
int L = N/M;
double xn[3] = {0, 0, 0};
double delta = 2.56; //step per metropolis
double delta2 = 6.95;
double delta3= 4.8;
double a0 = 1.;
double xn1[3] = {0, 0, 0};

double MODPSI(double x, double y, double z){ //funzione d'onda atomo idrogeno, nlm=100
	double u = pow(a0, -3)*exp(-2*sqrt(x*x+y*y+z*z)/a0)/M_PI;
	return u;
}

double MODPSI2(double x, double y, double z){ //nlm = 210
	double r = sqrt(x*x+y*y+z*z);
	double u = pow(a0, -5)*exp(-r/a0)*z*z/(32*M_PI);
	return u;
}

int main(){
	Random rnd;
	Blocchi  bl(N, M, "output.dat");
	Blocchi  bl2(N, M, "output2.dat");
	Blocchi bl3(N,M, "output3.dat");
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
	
	//metropolis step START	
	int nacc =0, nthrow= 0;
	for(int p=0; p < M; p++){
		xn[0]=0.;xn[1]=0.;xn[2]=0.;
		for(int q = 0; q<L; q++){
			int d = int(3*rnd.Rannyu());
			for(int i = 0; i < 3; i++){
				xn1[i] = xn[i];
			}
				xn1[d]=rnd.Rannyu(xn[d]-delta, xn[d]+delta);
				double a=min(1., MODPSI(xn1[0], xn1[1], xn1[2])/MODPSI(xn[0], xn[1], xn[2]));
				nthrow++;
				if(rnd.Rannyu() <= a){ //accept
					xn[d]=xn1[d];
					nacc++;
				}
				else { //reject
					xn1[d]=xn[d];
				} 
			bl.addThrow(sqrt(xn[0]*xn[0]+xn[1]*xn[1]+xn[2]*xn[2]));
		}
	}
	cout << "accettazione: " << (double)nacc/nthrow << endl;
	bl.saveoutput();
	nacc =0; nthrow= 0;
	for(int p=0; p < M; p++){
		xn[0]=0.;xn[1]=0.;xn[2]=0.;
		for(int q = 0; q<L; q++){
			int d = int(3*rnd.Rannyu());
			for(int i = 0; i < 3; i++){
				xn1[i] = xn[i];
			}
				xn1[d]=rnd.Rannyu(xn[d]-delta2, xn[d]+delta2);
				double a=min(1., MODPSI2(xn1[0], xn1[1], xn1[2])/MODPSI2(xn[0], xn[1], xn[2]));
				nthrow++;
				if(rnd.Rannyu() <= a){ //accept
					xn[d]=xn1[d];
					nacc++;
				}
				else { //reject
					xn1[d]=xn[d];
				}  
			bl2.addThrow(sqrt(xn[0]*xn[0]+xn[1]*xn[1]+xn[2]*xn[2]));
		}
	}
	cout << "accettazione: " << (double)nacc/nthrow << endl;
	bl2.saveoutput();
	//transizione gaussiana
	nacc =0; nthrow= 0;
	for(int p=0; p < M; p++){
		xn[0]=0.;xn[1]=0.;xn[2]=0.;
		for(int q = 0; q<L; q++){
			int d = int(3*rnd.Rannyu());
			for(int i = 0; i < 3; i++){
				xn1[i] = xn[i];
			}
				xn1[d]=rnd.Gauss(xn[d], delta3);
				double a=min(1., MODPSI2(xn1[0], xn1[1], xn1[2])/MODPSI2(xn[0], xn[1], xn[2]));
				nthrow++;
				if(rnd.Rannyu() <= a){ //accept
					xn[d]=xn1[d];
					nacc++;
				}
				else { //reject
					xn1[d]=xn[d];
				}  
			bl3.addThrow(sqrt(xn[0]*xn[0]+xn[1]*xn[1]+xn[2]*xn[2]));
		}
	}
	cout << "accettazione: " << (double)nacc/nthrow << endl;
	bl3.saveoutput();
	//test partenza lontanto dall'origine
	cout << "testing behaviour when starting far from origin: (0,100,0)" << endl;
	ofstream pos("pos.xyz");
	xn[0]=0.;xn[1]=100.;xn[2]=0.;
	for(int q = 0; q<40000; q++){
		int d = int(3*rnd.Rannyu());
		for(int i = 0; i < 3; i++){
			xn1[i] = xn[i];
		}
			xn1[d]=rnd.Rannyu(xn[d]-delta2, xn[d]+delta2);
			double a=min(1., MODPSI2(xn1[0], xn1[1], xn1[2])/MODPSI2(xn[0], xn[1], xn[2]));
			if(rnd.Rannyu() <= a){ //accept
				xn[d]=xn1[d];
			}
			else { //reject
				xn1[d]=xn[d];
			}  
			pos << xn[0] << "     " << xn[1] << "     " << xn[2] << endl;
	}
	pos.close();
  rnd.SaveSeed();
  return 0;
}


