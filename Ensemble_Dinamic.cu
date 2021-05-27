#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>

#include <sys/stat.h>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/openmp/openmp.hpp>
#include <boost/random.hpp>

typedef double prec;

#include "Ensemble_param.h"
#include "Ensemble_connections.h"
#include "Ensemble_Dinamic.h"


prec Dinamic::to_odeint::interaction_sum_all(int id, const std::vector<prec> &x)
{
	prec sum=0;
	for (int i = 0; i < N_t_odeint; ++i)
	{
		sum=sum+A_t_odeint(id,i)*sin(x[2*i]-x[2*id])/N_t_odeint;
	}
	return sum*K_t_odeint;
}

prec Dinamic::to_odeint::interaction_sum_chain(int id,const std::vector<prec> &x)
{
			prec sum=0;
		    if(id>0 && id<N_t_odeint-1)
		    {
		    	sum=sum+A_t_odeint(id,id+1)*sin(x[2*(id+1)]-x[2*id])/N_t_odeint;
		    	sum=sum+A_t_odeint(id,id-1)*sin(x[2*(id-1)]-x[2*id])/N_t_odeint;
		    	return sum*K_t_odeint;
 		    }
 		    if(id==0)
 		    {
		    	sum=sum+A_t_odeint(id,id+1)*sin(x[2*(id+1)]-x[2*id])/N_t_odeint;
		    	return sum*K_t_odeint;
		    }
		    if(id==N_t_odeint-1)
		    {
		    	sum=sum+A_t_odeint(id,id-1)*sin(x[2*(id-1)]-x[2*id])/N_t_odeint;
		    	return sum*K_t_odeint;
		    }  
			return sum*K_t_odeint;
}

prec Dinamic::to_odeint::interaction_sum(int id, const std::vector<prec> &x)
{
	if(A_t_odeint.get_type()=="chain")
	{
		return interaction_sum_chain(id,x);
	}
	if(A_t_odeint.get_type()=="global" || A_t_odeint.get_type()=="custom")
	{
		return interaction_sum_all(id,x);
	}
}


void Dinamic::print_params()
{
	std::ofstream txtOut;
	txtOut.open("params.txt");
	txtOut.precision(8);
	txtOut << N << " " << K << std::endl;
	txtOut.close();
	A.print();
   	I.print();
   	F.print();
   	G.print();
   	W.print();
}

void Dinamic::print_params_to_console()
{
	A.print_to_console();
	I.print_to_console();
	F.print_to_console();
	G.print_to_console();
	W.print_to_console();
}

void Dinamic::generate()
{
	A.generate(A.get_type());
   	I.generate();
   	F.generate();
   	G.generate();
   	W.generate();
}