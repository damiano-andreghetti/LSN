#include "GA.h"
#define _USE_MATH_DEFINES 
using namespace std;

double distance(double x1, double y1, double x2, double y2){
	double d = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
	return d;
}

int GA::GenPBC(int i){
	int r=i;
	if(i>Ng-1) r = i%Ng;
	else if(i < 0) r = Ng-(abs(i)%Ng);
}
GA::GA(int N1, int N2){
	Np = N1;
	Ng = N2;
	Pop.resize(Np);
	for(int i = 0; i < Np; i++) Pop[i].resize(Ng);
	cities.resize(Ng);
	for(int i = 0; i < Ng; i++) cities[i].resize(2);
	Length.resize(Np);
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
}

GA::~GA(){}

void GA::GenerateCities(int N, int mode){
	if(mode==0){
		for(int i = 0; i < N; i++){
			double teta = rnd.Rannyu(0., 2.*M_PI);
			double xc = cos(teta);
			double yc = sin(teta);
			//cout << xc << "," << yc << endl;
			//cities[i] = vector<double>(2);
			cities[i][0]=xc;
			cities[i][1]=yc;
		}
	}
	else if(mode==1){
		for(int i = 0; i < N; i++){
			cities[i][0]=rnd.Rannyu(-1., 1.);
			cities[i][1]=rnd.Rannyu(-1., 1.);
		}
	}
	rnd.SaveSeed();
}

void GA::SaveBest(int n){
	string fname="best_individual";
	fname += to_string(n);
	fname += ".dat";
	ofstream best(fname);
	for(int i = 0; i < Ng; i++){ // export best individual
		best << Pop[0][i] << endl; 
	}
	best << "cities: " <<endl;
	for(int i=0; i< cities.size(); i++){ //export cities positions
		best << cities[i][0] << " " << cities[i][1] << endl;
	}
	best.close();
}

void GA::Stats(){
	stats.push_back(Length[0]);
	double av = 0;
	for(int i = 0; i < Np/2; i++){
		av+=Length[i];
	}
	av= av/(Np/2);
	avestats.push_back(av);
}

void GA::SaveStats(int n){
	string fname="best_stats";
	fname += to_string(n);
	fname += ".dat";
	ofstream fstats(fname);
	for(int i = 0; i <stats.size(); i++){
		fstats << stats[i] << endl;
	}
	fstats.close();
	string fname2="average_stats";
	fname2 += to_string(n);
	fname2 += ".dat";
	ofstream fave(fname2);
	for(int i = 0; i < avestats.size(); i++){
		fave << avestats[i] << endl;
	}
	fave.close();
}
void GA::GeneratePop(){
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

void GA::EvaluatePop(){
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

void GA::SortPop(){
	bool sorted;
	do{
		sorted=1;
		for(int i =1; i < Np; i++){
			if(Length[i]<Length[i-1]){ //swap elements
				sorted=0;
				double v = Length[i];
				Length[i]=Length[i-1];
				Length[i-1]=v;
				Pop[i].swap(Pop[i-1]	);
			}
		}
	}
	while(!sorted);
	cout << endl << "cost of the best individual of this generation:" << endl;
	cout << Length[0] << endl;
}

int GA::Selection(){
	return int(Np*pow(rnd.Rannyu(),2.));
}

void GA::Mutate(int m, int i){ // m = type of mutation. i = individual to mutate
	if(m==0){ //permutate 2 genes
		int a = int(Ng*rnd.Rannyu());
		int b = GenPBC(int(Ng*rnd.Rannyu()));
		int t =	Pop[i][a];
		Pop[i][a] = Pop[i][b];
		Pop[i][b] = t;
	}
	else if(m==1){ //shift genes
		int s = int(Ng*0.5*rnd.Rannyu());
		rotate(Pop[i].begin(), Pop[i].begin()+s, Pop[i].end());			
	}
	else if(m==2){ // shift a group of genes
		int s = 1+int(Ng*0.4*rnd.Rannyu()); //shift parameter
		int d = 1+int(Ng*0.2*rnd.Rannyu()); //dimension of the group to shift
		int b = int(Ng*rnd.Rannyu()); //beginning index of the group to shift
		int bk[d];
		for(int k = 0; k < d; k++){
			bk[k]=Pop[i][GenPBC(b+k)];
		}
		for(int k = 0; k<s; k++){
			Pop[i][GenPBC(b+k)]=Pop[i][GenPBC(b+d+k)];
		}
		for(int k=0; k < d; k++){
			Pop[i][GenPBC(b+s+k)]=bk[k];
		}
	}
	else if(m==3){ //permutation of a group of genes
		int a = int(Ng*rnd.Rannyu());
		int d = 1+int(0.2*Ng*rnd.Rannyu());
		int b = GenPBC(a+ d + int(Ng*0.2*rnd.Rannyu()));	
		int buff;
		for(int k = 0; k < d; k++){
			buff=Pop[i][GenPBC(a+k)];
			Pop[i][GenPBC(a+k)]=Pop[i][GenPBC(b+k)];
			Pop[i][GenPBC(b+k)]=buff;
		}	
	}
	else if(m==4){ //inversion of a group of genes
		int a = int(Ng*rnd.Rannyu());
		int d = 1+int(0.2*Ng*rnd.Rannyu());
		int bk[d];
		for(int k = 0; k < d; k++){
			bk[k]=Pop[i][GenPBC(a+k)];
		}
		for(int k = 0; k < d; k++){
			Pop[i][GenPBC(a+d-k)]=bk[k];
		}
	}
	else{
	cout << "incorrect mutation index!!" << endl;
	}
}

void GA::PopMutation(double p){ //probability to get a mutation;
	for(int i=0; i < Np; i++){
		for(int m = 0; m < Nm; m++){ //cycle over all mutation
			if(rnd.Rannyu()<p) Mutate(m, i);
		}
	}
}
bool GA::IsInVector(int a, vector<int> vec){
	bool found = false;
	for(int i = 0; i < vec.size(); i++){
		if(vec[i]==a) found=true;
	}
	return found;
}
void GA::NewGeneration(){ //Crossover function
	vector<vector<int> > Offspring(Np);
	for(int i = 0; i < Np; i++) Offspring[i].resize(Ng);
	for(int j = 0; j< Np; j+=2){
		//select two parents 
		int a = Selection();
		int b = Selection();
		int cut = Ng*rnd.Rannyu();
		for(int k = 0; k < cut; k++){
			Offspring[j][k]=Pop[a][k];
			Offspring[j+1][k]=Pop[b][k];
		}
		int o1 = cut, o2 = cut;
		for(int k = 0; k < Ng; k++){
			if(!IsInVector(Pop[b][k],Offspring[j])){
				Offspring[j][o1]=Pop[b][k];
				o1++;
			}
		}
		for(int k = 0; k < Ng; k++){
			if(!IsInVector(Pop[a][k],Offspring[j+1])){
				Offspring[j+1][o2]=Pop[a][k];
				o2++;
			}
		}
	}
	Pop.swap(Offspring);
}
