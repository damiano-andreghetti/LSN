#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "blocchi.h"

#define _USE_MATH_DEFINES 
using namespace std;


unsigned long int N = 2000000;
unsigned long int M = 100;
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
	//double u = exp(-(x-m)*(x-m)/(2*s*s))+exp(-(x+m)*(x+m)/(2*s*s));
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
/*
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
*/
	double mu =0.8;
	ofstream samp("samples.dat");
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
			samp << xn << endl;
			double stima_H = potential(xn);
			bl.addThrow(stima_H);
			}
	}
	cout << acc*100./N << endl;
	samp.close();
  rnd.SaveSeed();
  return 0;
}


