#ifndef __Blocchi__
#define __Blocchi__

#include <cstdlib>
#include <iostream>
#include <fstream>
class Blocchi {

private:
  int m_N, m_M, m_L, m_i, m_j; //numero punti, numero blocchi, lunghezza blocchi, indice blocco, indice punto nel blocco 
	double *av, *av2;
	char* m_filename;
protected:

public:
  // constructors
	Blocchi(int, int, char*);
  // destructor
  ~Blocchi();
  // methods
  void  addThrow(double);
	void saveoutput();
};

#endif // __Blocchi__
