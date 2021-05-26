#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>

#include <sys/stat.h>
#include <boost/numeric/odeint.hpp>
#include <boost/random.hpp>

typedef double prec;

#include "Ensemble_param.h"
#include "Ensemble_connections.h"
#include "Ensemble_Dinamic.h"


int main()
{
 	boost::mt19937 rng(static_cast<unsigned int>(std::time(0)));

	Dinamic P(9 ,1 ,rng ,1 ,0.1 ,2 ,0.2 ,3 ,0.3 ,4 ,0.3,"global");
	//P.print_params();
	/*Dinamic P(rng);*/
	P.print_params_to_console();
	P.generate();
	P.print_params_to_console();

	return 0;
}