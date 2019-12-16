#include "blocchi.h"
#include <cmath>

using namespace std;

Blocchi::Blocchi(int N, int M, char* filename){
	m_N=N;
	m_M=M;
	m_L=m_N/m_M;
	m_i=0;
	m_j=0;
	m_filename=filename;
	av = new double [m_M];
	av2 = new double [m_M];
}
Blocchi::~Blocchi(){}

void Blocchi::addThrow(double u){
	if(m_j==m_L-1){av[m_i]/=m_L; av2[m_i]=av[m_i]*av[m_i]; m_i++; m_j=0;}
	if(m_i==m_M){
		saveoutput();	
		return;
	}
	if(m_j==0) {av[m_i]=0.; av2[m_i]=0.;}
	av[m_i]+=u;
	m_j++;
}

void Blocchi::saveoutput(){
	ofstream output;
	output.open(m_filename);
	for(int  i = 0; i < m_M; i++){
		double s=0., sq=0., e=0.;
		for(int j = 0; j<i+1; j++){
			s+=av[j];
			sq+=av2[j];
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
}
