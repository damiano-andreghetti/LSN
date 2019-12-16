#ifndef __Blocchi__
#define __Blocchi__

#include <cstdlib>
#include <iostream>
#include <fstream>
class Blocchi {

private:
  unsigned long int m_N, m_M, m_L, m_i, m_j; 
	bool m_full;
//numero punti, numero blocchi, lunghezza blocchi, indice blocco, indice punto nel blocco 
	double *av, *av2;
	char* m_filename;
protected:

public:
  // constructors
	Blocchi(unsigned long int, unsigned long int, char*);
  // destructor
  ~Blocchi();
  // methods
  void  addThrow(double);
	void saveoutput();
};

#endif // __Blocchi__
