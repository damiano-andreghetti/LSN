/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>
#include <algorithm> 

#define _USE_MATH_DEFINES 

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

	ofstream output("output.dat");
	ofstream output2("output2.dat");
	ofstream output3("output3.dat");
	ofstream output4("output4.dat");

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
	int N = 1000000; //punti totale
	int M = 100; //numero blocchi
	int L = N/M; //lunghezza singolo blocco
	float* av = new float[M];
	float* avq = new float[M];
	//es 3.1 Black-Scholes theory
	float T=1.; //delivery time
	float S0=100.; //asset price at t=0 
	float K=100.;	//strike-price
	float r=0.1;	//risk-free interest rate
	float sigma=0.25;	//volatility
	//codice per esportare risultato in funzione dei blocchi e con errore 
	//call option direct s(T) sampling
	cout << "direct sampling S(T) for call option" << endl;
	for(int i=0; i<M;i++){ //ciclo blocchi
		av[i]=0.;
		avq[i]=0.;
		for(int j=0; j<L; j++){ 
			double z = rnd.Gauss(0,1);
			double S = S0*exp((r-sigma*sigma/2)*T+sigma*z*sqrt(T));
			av[i]+=exp(-r*T)*max(0., S-K);
		}	
		av[i]/= L;
		avq[i] = av[i]*av[i];
	}
	for(int i=0; i<M; i++){
		float s=0.;
		float sq=0.;
		float e=0.;
		for(int j=0; j<i+1; j++){
			s+=av[j];
			sq+=avq[j];
		}
		s/=(i+1);
		sq/=(i+1);
		if(i==0) e=0;
		else{
			e=sqrt((abs(sq-s*s))/i);
		}
		output << i << " " << s << " " << e << endl;
	}
	output.close();
	cout << "direct sampling S(T) for put option" << endl;
	for(int i=0; i<M;i++){ //ciclo blocchi
		av[i]=0.;
		avq[i]=0.;
		for(int j=0; j<L; j++){ 
			double z = rnd.Gauss(0,1);
			double S = S0*exp((r-sigma*sigma/2)*T+sigma*z*sqrt(T));
			av[i]+=exp(-r*T)*max(0., K-S);
		}	
		av[i]/= L;
		avq[i] = av[i]*av[i];
	}
	for(int i=0; i<M; i++){
		float s=0.;
		float sq=0.;
		float e=0.;
		for(int j=0; j<i+1; j++){
			s+=av[j];
			sq+=avq[j];
		}
		s/=(i+1);
		sq/=(i+1);
		if(i==0) e=0;
		else{
			e=sqrt((abs(sq-s*s))/i);
		}
		output2 << i << " " << s << " " << e << endl;
	}
	output2.close();
	cout << "discrete sampling S(T) through 100 intervals for call option" << endl;
	for(int i=0; i<M;i++){ //ciclo blocchi
		if (i%10==0) cout << "block n°" << i << " of " << M <<endl;
		av[i]=0.;
		avq[i]=0.;
		for(int j=0; j<L; j++){ 
			double S = S0; 
			for(int k = 0; k < 100; k++){
				double z = rnd.Gauss(0,1);
				S=S*exp((r-sigma*sigma/2)*(T/100.)+sigma*z*sqrt(T/100.));
			}
			av[i]+=exp(-r*T)*max(0., S-K);
		}	
		av[i]/= L;
		avq[i] = av[i]*av[i];
	}
	for(int i=0; i<M; i++){
		float s=0.;
		float sq=0.;
		float e=0.;
		for(int j=0; j<i+1; j++){
			s+=av[j];
			sq+=avq[j];
		}
		s/=(i+1);
		sq/=(i+1);
		if(i==0) e=0;
		else{
			e=sqrt((abs(sq-s*s))/i);
		}
		output3 << i << " " << s << " " << e << endl;
	}
	output3.close();
	cout << "discrete sampling S(T) through 100 intervals for put option" << endl;
	for(int i=0; i<M;i++){ //ciclo blocchi
		if (i%10==0) cout << "block n°" << i << " of " << M <<endl;
		av[i]=0.;
		avq[i]=0.;
		for(int j=0; j<L; j++){ 
			double S = S0; 
			for(int k = 0; k < 100; k++){
				double z = rnd.Gauss(0,1);
				S=S*exp((r-sigma*sigma/2)*(T/100.)+sigma*z*sqrt(T/100.));
			}
			av[i]+=exp(-r*T)*max(0., K-S);
		}	
		av[i]/= L;
		avq[i] = av[i]*av[i];
	}
	for(int i=0; i<M; i++){
		float s=0.;
		float sq=0.;
		float e=0.;
		for(int j=0; j<i+1; j++){
			s+=av[j];
			sq+=avq[j];
		}
		s/=(i+1);
		sq/=(i+1);
		if(i==0) e=0;
		else{
			e=sqrt((abs(sq-s*s))/i);
		}
		output4 << i << " " << s << " " << e << endl;
	}
	output4.close();
  rnd.SaveSeed();
  return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
