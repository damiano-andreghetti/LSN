/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
const int m_props=1000;
double walker[m_props];
int n_props;
int indexave = 0;
int iv,ik,it,ie, igofr;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_pres, bin_size, nbins;
double *av_pot, *av2_pot, *av_kin, *av2_kin, *av_tot, *av2_tot, *av_temp, *av2_temp, *av_pres, *av2_pres;
int cb = 0; //current block counter;

// averages
double acc,att;
double glob_av[m_props],glob_av2[m_props];
double blk_av[m_props],blk_norm;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed, newstart, isaveconf, nblocks, L, imeasure;
double delta;

//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Gdirave(void);
void Accumualte(void);
void Reset(int);
void Measure(int);
double Force(int, int);
double Pbc(double);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
