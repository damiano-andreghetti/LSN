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

using namespace std;

double PBC(double x){ //periodic boundary condiions for buffon's experiment simulation
	double r = x;
	if(x>1.) r = x-1.;
	else if(x<0.) r = 1.+x;
	return r;
}
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
	ofstream outputX("outputX.dat");

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
	//codice per esportare risultato in funzione dei blocchi e con errore	
	for(int i=0; i<M;i++){ //ciclo blocchi
		av[i]=0.;
		avq[i]=0.;
		for(int j=0; j<L; j++){ 
			float u = rnd.Rannyu();
			av[i]+=u;
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
	//test del chi quadro, per la verifica delle ipotesi
	int MX= 100;
	int nX= 10000;
	int LX = nX/MX;
	double d = 1./MX;
	for(int j = 0; j < 100; j++){ // eseguo il test 100 volte
		double X = 0.;
		for(int i=0; i<MX;i++){ // divido [0,1] in 100 sottointervalli
			double nk = 0;
			for(int k = 0; k < LX; k++){
				nk += rnd.Rannyu(i*d, (i+1)*d);
			}
			nk=nk/LX;
			X += (nk-nX/MX)*(nk-nX/MX)/(nX/MX);
		}
		X=X/MX;
		outputX << j << " " << X << endl;
	} 
	outputX.close();
	//es1.2
	for(int i=0; i<M;i++){ //ciclo blocchi
		av[i]=0.;
		avq[i]=0.;
		for(int j=0; j<L; j++){ 
			float u = rnd.Rannyu();
			av[i]+=(u-0.5)*(u-0.5);
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
	//es 1.2
	for(int i=0; i<10000; i++){ //test distribuzioni
		float x1 = rnd.Exponential(1.);
		float x2 = rnd.CauchyLorentz(0.,1.);
		output3 << x1 << " " << x2 << endl;
	}
	output3.close();
	//es 1.2 vero e proprio
	ofstream histo("histo.dat");
	for(int i = 0; i < 10000; i++){
		double s[4] = {0}; //S1, S2, S10, S100 for standard dice
		double e[4] = {0}; // " "  for exp. dice
		double l[4] = {0}; // " "  for Lorentz dice
		s[0] += rnd.Rannyu();
		e[0] += rnd.Exponential(1.);
		l[0] += rnd.CauchyLorentz(0., 1.);
		s[1] += rnd.Rannyu() + rnd.Rannyu(); s[1]/=2;
		e[1] += rnd.Exponential(1.) + rnd.Exponential(1.); e[1]/=2;
		l[1] += rnd.CauchyLorentz(0., 1.); l[1]/=2;
		for(int k = 0; k < 10; k++){
			s[2]+=rnd.Rannyu();
			e[2]+=rnd.Exponential(1.);
			l[2]+=rnd.CauchyLorentz(0.,1.);
		}
		s[2]/=10; e[2]/=10; l[2]/=10;
		for(int k = 0; k < 100; k++){
			s[3]+=rnd.Rannyu();
			e[3]+=rnd.Exponential(1.);
			l[3]+=rnd.CauchyLorentz(0.,1.);
		}
		s[3]/=100; e[3]/=100; l[3]/=100;
		histo << s[0] << " " << s[1] << " " << s[2] << " " << s[3] << " " ;
		histo << e[0] << " " << e[1] << " " << e[2] << " " << e[3] << " " ;
		histo << l[0] << " " << l[1] << " " << l[2] << " " << l[3] << endl;
	}
	histo.close();	
	//es 1.3
	ofstream buff("outputbuff.dat");
	double dist = 1.;
	double Ln = 0.782;	
	L=1000;
	for(int i=0; i<M;i++){ //ciclo blocchi
		av[i]=0.;
		avq[i]=0.;
		int Nthr=0, Nhit=1;
		for(int j=0; j<L; j++){ 
			double x0,dy, x1;
			Nthr++;
			x0 = rnd.Rannyu()*dist;
			dy = rnd.Rannyu(-Ln, Ln);
			x1 = sqrt((Ln*Ln)-(dy*dy));
			double u = rnd.Rannyu();
			if(u<0.5) x1+=x0;
			else x1-=x0;
			if(x1<0 ||x1>dist) Nhit++;
			cout << 2*Ln*Nthr/(Nhit*dist)<< endl;
			av[i]+= 2*Ln*Nthr/(Nhit*dist);
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
		buff << i << " " << s << " " << e << endl;
	}
	buff.close();		
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
****************************************************************/
