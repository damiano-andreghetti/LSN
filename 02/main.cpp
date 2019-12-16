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
	int N = 100000; //punti totale
	int M = 100; //numero blocchi
	int L = N/M; //lunghezza singolo blocco
	float* av = new float[M];
	float* avq = new float[M];
	//es 2.1 uniforme
	//codice per esportare risultato in funzione dei blocchi e con errore
	for(int i=0; i<M;i++){ //ciclo blocchi
		av[i]=0.;
		avq[i]=0.;
		for(int j=0; j<L; j++){ 
			float u = rnd.Rannyu();
			av[i]+=(0.5*M_PI)*cos(M_PI*u*0.5);
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
	//es2.1 importance sampling d(x)=-2x+2
	for(int i=0; i<M;i++){ //ciclo blocchi
		av[i]=0.;
		avq[i]=0.;
		for(int j=0; j<L; j++){ 
			float u = rnd.samplear();
			av[i]+=M_PI*cos(M_PI*u*0.5)/(-4*u+4);
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
	//es 2.2 random walk on 3D lattice (discrete)
	//100000 points, plot sqrt(|r|^2) in function of number of steps (100 steps for each block)
	ofstream disc("randomwalkdisc.dat");
	int T = 10000;
	int steps = 100;
	double media[steps] = {0};
	double mediaq[steps] = {0};
	double a = 1.;
	for(int i = 0; i < T; i++){ //repeat the simulation 10000 times 
		double pos[3] = {0,0,0};
		double r = 0;
		for(int j = 0; j < steps; j++){ //100 steps for each random walk
			int d = int(3*rnd.Rannyu()); //select random direction between x,y,z
			int ad = 1;
			double u = rnd.Rannyu();
			if(u<0.5) ad = -1;
			pos[d]+=ad*a;
			r = sqrt(abs(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]));
			media[j]+=r;
			mediaq[j]+=r*r;
		}
	}
	double s =0.,sq = 0., e = 0.; 
	for(int j = 0; j < steps; j++){
		media[j]/=T;
		mediaq[j]/=T;
		s = ((s*j) + media[j])/(j+1.);
		sq = ((sq*j) + media[j])/(j+1.);
		if (j == 0) e = 0.;
		else e=sqrt((abs(sq-s*s))/(j));
		disc << j << " " << s << " " << e << endl;
	}
	disc.close();
	//random walk continua
	ofstream cont("randomwalkcont.dat");
	for(int j = 0; j < steps; j++){
		media[j]=0.;
		mediaq[j]=0.;
	}
	for(int i = 0; i < T; i++){ //repeat the simulation 10000 times 
		double pos[3] = {0};
		double r = 0;
		for(int j = 0; j < steps; j++){ //100 steps for each random walk
			double theta = rnd.Rannyu(0., M_PI);
			double phi = rnd.Rannyu(0., 2.*M_PI);
			pos[0]+=a*sin(theta)*cos(phi);
			pos[1]+=a*sin(theta)*sin(phi);
			pos[2]+=a*cos(theta);
			r = sqrt(abs(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]));
			media[j]+=r;
			mediaq[j]+=r*r;
		}
	}
	s =0.;sq = 0.; e = 0.; 
	for(int j = 0; j < steps; j++){
		media[j]/=T;
		mediaq[j]/=T;
		s = ((s*j) + media[j])/(j+1.);
		sq = ((sq*j) + media[j])/(j+1.);
		if (j == 0) e = 0.;
		else e=sqrt((abs(sq-s*s))/(j));
		cont << j << " " << s << " " << e << endl;
	}
	cont.close();
	/*
	ofstream test("test.dat");
	for(int i =0; i<10000; i++) test << rnd.samplear()<< endl;
	test.close();
	*/
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
