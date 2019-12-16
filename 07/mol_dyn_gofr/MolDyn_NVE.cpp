/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <iomanip>
#include "MolDyn_NVE.h"
#include "random.h"


using namespace std;

int main(){ 
  Input();             //Inizialization
  int nconf = 1;
	int mstep = 0;
  for(int istep=1; istep <= nstep; ++istep){
     Move();           //Move particles with Verlet algorithm
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%isaveconf == 0){
        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
				mstep++;
     }
		 if(istep%imeasure==0){
		 Measure(istep/imeasure);
		 }
  }
	Gdirave();
  ConfFinal();         //Write final configuration to restart

  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf,ReadOldConf;
  double ep, ek, pr, et, vir;
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

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
	ReadInput >> newstart;
	ReadInput >> isaveconf;
	ReadInput >> nblocks;
	ReadInput >> imeasure;

	L = nstep/(nblocks*imeasure); //block length
	av_pot = new double[nblocks];
	av2_pot = new double[nblocks];
	av_kin = new double[nblocks];
	av2_kin = new double[nblocks];
	av_tot = new double[nblocks];
	av2_tot = new double[nblocks];
	av_temp = new double[nblocks];
	av2_temp = new double[nblocks];
	av_pres = new double[nblocks];
	av2_pres = new double[nblocks];
  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  cout << "Number of blocks = " << nblocks << " of " << L << " steps each" << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables
 //measurement of g(r)
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;

	if(newstart == 1){
	//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; i++){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();
//Prepare initial velocities
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; i++){
     vx[i] = rnd.Rannyu() - 0.5;
     vy[i] = rnd.Rannyu() - 0.5;
     vz[i] = rnd.Rannyu() - 0.5;
     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; i++){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;
   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   for (int i=0; i<npart; i++){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = x[i] - vx[i] * delta;
     yold[i] = y[i] - vy[i] * delta;
     zold[i] = z[i] - vz[i] * delta;
   }
	}
/*
Se invece possiedo un file delle ultime due configurazioni posso calcolarmi le velocita 
di coseguenza, riscalandole poi per adattarle alla temperatura
*/
	if(newstart!=1){
//Read initial configuration
  cout << "Read initial configuration from file config.final " << endl << endl;
  ReadConf.open("config.final");
  for (int i=0; i<npart; i++){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();
		ReadOldConf.open("oldconfig.final");
	  for (int i=0; i<npart; i++){
    	ReadOldConf >> xold[i] >> yold[i] >> zold[i];
   		xold[i] = xold[i] * box;
    	yold[i] = yold[i] * box;
    	zold[i] = zold[i] * box;
  	}
  	ReadConf.close(); 
		//faccio fare un passo al sistema, trovo r(t+dt)
		//lo uso per stimare le v(t), le riscalo e ricalcolo r(t-dt) di conseguenza
		double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];
		double sumv2 = 0., fs;
		for(int i=0; i<npart; i++){ //Force acting on particle i
	    fx[i] = Force(i,0);
	    fy[i] = Force(i,1); 	 		 	
	    fz[i] = Force(i,2);
	  }
		for(int i =0; i<npart; i++){
		xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );
		vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);
		sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
		x[i]=xnew;
		y[i]=ynew;
		z[i]=znew;
		}
  	sumv2 /= (double)npart;
		sumv2=sqrt(sumv2);
		cout << "velocità media particelle precedenti: " << sumv2 <<endl;
		cout << "velocità media attesa con equipartizione: " << sqrt(3.*temp) << endl;
		fs = sqrt(3 * temp)/sumv2;   // fs = velocity scale factor 
		cout << "fattore per riscalare le vocita': " << fs << endl;
		for (int i=0; i<npart; i++){
	    vx[i] *= fs;
	    vy[i] *= fs;
	    vz[i] *= fs;
			sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
	    xold[i] = Pbc(x[i] - 2. * vx[i] * delta);
	    yold[i] = Pbc(y[i] - 2. * vy[i] * delta);
	    zold[i] = Pbc(z[i] - 2. * vz[i] * delta);
   }
		sumv2/= (double)npart;
		sumv2=sqrt(sumv2);
		cout << "velocità media particelle precedenti riscalate secondo equipartizione" << sumv2 << endl;
ReadConf.open("config.final");
  for (int i=0; i<npart; i++){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();
	}
  rnd.SaveSeed();
	return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}
void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
}


void Measure(int mstep){ //Properties measurement
  int bin;
  double v, t, vij, p, r, gdir;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Gofr;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;
	p = 0.0;	
//reset the hystogram of g(r)
  for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;
//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);
	bin = int(dr/bin_size);
		walker[igofr+bin]++;
     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
			p += (pow(dr, -12.) - 0.5*pow(dr, -6.));
//Potential energy
       v += vij;
     }
    }    
	Accumulate();      
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
    stima_pot = v/(double)npart; //Potential energy
    stima_kin = t/(double)npart; //Kinetic energy
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total enery
 		stima_pres = rho*stima_temp + 48*p*(1/(3.*vol));
		//metodo blocchi 
		av_pot[cb]+=stima_pot;
		av_kin[cb]+=stima_kin;
		av_tot[cb]+=stima_etot;
		av_temp[cb]+=stima_temp;
		av_pres[cb]+=stima_pres;
		if(mstep%L==0){
			av_pot[cb] /=(double)L;
			av_kin[cb] /=(double)L;
			av_tot[cb] /=(double)L;
			av_temp[cb] /=(double)L;
			av_pres[cb] /=(double)L;
			av2_pot[cb] =av_pot[cb]*av_pot[cb];
			av2_kin[cb] =av_kin[cb]*av_kin[cb];
			av2_tot[cb] =av_tot[cb]*av_tot[cb];
			av2_temp[cb] =av_temp[cb]*av_temp[cb];
			av2_pres[cb] =av_pres[cb]*av_pres[cb];	
			for(int i =0; i<nbins; i++){
				double norm =0.;
				r=i*bin_size;
				norm=rho*npart*(4*M_PI*(pow(r+bin_size, 3.)-pow(r, 3.)))/3;
				gdir = blk_av[igofr+i]/blk_norm/norm;
				glob_av[igofr+i]+=gdir;
				glob_av2[igofr+i]+=gdir*gdir;
			}
			Reset(cb);
			cb++;
		}
    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    return;
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;
	
  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
	ofstream WriteOldConf;

  cout << "Print configuration before the final one to file oldconfig.final " << endl << endl;
  WriteOldConf.open("oldconfig.final");

  for (int i=0; i<npart; ++i){
    WriteOldConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteOldConf.close();
	string s;
	ifstream in;
	in.open("ave_epot.dat");
	if(in.is_open()){
	while(!in.eof()) {
		getline(in, s);
		indexave ++;	
	}
	cout << "appending file containing averages to previous simulation" << endl;	
	}
	else { cout << "no previous simulation found" << endl;}
	in.close();
	indexave--;
	//analisi dati e creazione file con valori medi ed incertezze
	ofstream ave_epot, ave_ekin, ave_etot, ave_temp, ave_pres;
	ave_epot.open("ave_epot.dat",ios::app);
	ave_ekin.open("ave_ekin.dat",ios::app);
	ave_etot.open("ave_etot.dat",ios::app);
	ave_temp.open("ave_temp.dat",ios::app);
	ave_pres.open("ave_pres.dat",ios::app);
	for(int i=0; i<nblocks; i++){
		double s1=0., s2=0., s3=0., s4=0., s5=0.;
		double sq1=0., sq2=0., sq3=0., sq4=0., sq5=0.;
		double e1=0., e2=0., e3=0., e4=0., e5=0.;
		for(int j=0; j<i+1; j++){
			s5+=av_pres[j];
			sq5+=av2_pot[j];
			s1+=av_pot[j];
			sq1+=av2_pot[j];
			s2+=av_kin[j];
			sq2+=av2_kin[j];
			s3+=av_tot[j];
			sq3+=av2_tot[j];
			s4+=av_temp[j];
			sq4+=av2_temp[j];
		}
		s1/=(i+1);
		sq1/=(i+1);
		s2/=(i+1);
		sq2/=(i+1);
		s3/=(i+1);
		sq3/=(i+1);
		s4/=(i+1);
		sq4/=(i+1);
		s5/=(i+1);
		sq5/=(i+1);
		if(i==0) e1=0, e2=0, e3=0, e4=0, e5=0;
		else{
			e1=sqrt((abs(sq1-s1*s1))/i);
			e2=sqrt((abs(sq2-s2*s2))/i);
			e3=sqrt((abs(sq3-s3*s3))/i);
			e4=sqrt((abs(sq4-s4*s4))/i);
			e5=sqrt((abs(sq5-s5*s5))/i);
		}
		ave_epot << indexave+i << " " << s1 << " " << e1 << endl;
		ave_ekin << indexave+i << " " << s2 << " " << e2 << endl;
		ave_etot << indexave+i << " " << s3 << " " << e3 << endl;
		ave_temp << indexave+i << " " << s4 << " " << e4 << endl;
		ave_pres << indexave+i << " " << s5 << " " << e5 << endl;
	}
	ave_epot.close();
	ave_ekin.close();
	ave_etot.close();
	ave_temp.close();
	ave_pres.close();
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

void Gdirave(void){
	ofstream Gave;
	Gave.open("output.gave.0",ios::app);
	double err=0., r=0.;
	const int wd=12;
	cout << "exporting g(r) average with uncertainties" << endl;
	for(int i = 0; i<nbins; i++){
		r=i*bin_size;
		err=Error(glob_av[igofr+i], glob_av2[igofr+i], nblocks);
		Gave << setw(wd) << r << setw(wd) << glob_av[igofr+i]/(double)nblocks << setw(wd) << err << endl;
	}
	Gave.close();
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
