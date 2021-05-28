#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>
#include <limits>


#include <sys/stat.h>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/openmp/openmp.hpp>
#include <boost/random.hpp>

typedef double prec;

#include "Ensemble_param.h"
#include "Ensemble_connections.h"
#include "Ensemble_Dinamic.h"
#include "Ensemble_evolve.h"

void check_stationary(Evolve &e, Dinamic &P)
{
	bool statement=false;
	int count=0;
	while(statement==false)
	{
		e.clean();
		e.run(P,count*10000,(count+1)*10000,100);

		for (int i = 1; i < e.size(); ++i)
		{	
			prec convergence=e.convergence(i);
			if(convergence>0)
			{
				std::cout << i <<": " << convergence << std::endl;
			}
		}

		statement=e.did_converge(49);
		count++;

		e.reset_init();		
	}
	e.print();
	e.print_init();
}


int main()
{
 	boost::mt19937 rng(static_cast<unsigned int>(std::time(0)));

 	int N;
 	std::cout << "N: ";
 	std::cin >> N;

    int chunk_size = N/omp_get_max_threads();
    omp_set_schedule( omp_sched_static , chunk_size );

	Dinamic P(N ,20 ,rng ,1 ,0 ,5000 ,0 ,0.5 ,0 ,1 ,0,"chain"); //N, K, rng, I, sigm_I, F, sigm_F, G, sigm_G, W, sigm_W
	P.print_params();

	
	Evolve e(N, rng);
	//check_stationary(e, P);
	e.run(P,0,190000,100);
	e.print();






	return 0;
}