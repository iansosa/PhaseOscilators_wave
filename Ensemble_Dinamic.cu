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

typedef long double prec;

#include "Ensemble_param.h"
#include "Ensemble_connections.h"
#include "Ensemble_Dinamic.h"

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

prec Dinamic::to_odeint::h_interaction_sum_chain(int id,const std::vector<prec> &x)
{
			prec sum=0;
		    if(id>0 && id<N_t_odeint-1)
		    {

		    	sum=sum+A_t_odeint(id,id+1)*(x[2*(id+1)]-x[2*id])/N_t_odeint;
		    	sum=sum+A_t_odeint(id,id-1)*(x[2*(id-1)]-x[2*id])/N_t_odeint;		    	
		    	return sum*K_t_odeint;
 		    }
 		    if(id==0)
 		    {
		    	sum=sum+A_t_odeint(id,id+1)*(x[2*(id+1)]-x[2*id])/N_t_odeint;
		    	return sum*K_t_odeint;
		    }
		    if(id==N_t_odeint-1)
		    {
		    	sum=sum+A_t_odeint(id,id-1)*(x[2*(id-1)]-x[2*id])/N_t_odeint;
		    	return sum*K_t_odeint;
		    }  
			return sum*K_t_odeint;
}

prec Dinamic::to_odeint::interaction_sum(int id, const std::vector<prec> &x)
{
	switch(A_t_odeint.type_id)
	{
		case 0:
			return interaction_sum_all(id,x);
			break;
		case 1:
			return interaction_sum_chain(id,x);
			break;
		case 2:
			return interaction_sum_all(id,x);
			break;
		case 3:
			return h_interaction_sum_chain(id,x);
			break;
		case 4:
			return h_interaction_sum_chain(id,x);
			break;
		case 5:
			return h_interaction_sum_chain(id,x);
			break;
		case 6:
			return h_interaction_sum_chain(id,x);
			break;
		case 7:
			return 0;
			break;
		default :
			std::cout << "undefined ensemble type" << std::endl;
			return 0;
			break;
	}
}

prec Dinamic::to_odeint::force_sum(prec sum, int id, const std::vector<prec> &x,const prec t)
{
	switch(A_t_odeint.type_id)
	{
		case 0:
			return p_force_sum(sum,id,x,t);
			break;
		case 1:
			return p_force_sum(sum,id,x,t);
			break;
		case 2:
			return p_force_sum(sum,id,x,t);
			break;
		case 3:
			return h_force_sum(sum,id,x,t);
			break;
		case 4:
			return t_force_sum(sum,id,x,t);
			break;
		case 5:
			return s_force_sum(sum,id,x,t);
			break;
		case 6:
			return c_a_force_sum(sum,id,x,t);
			break;
		case 7:
			return solid_force_sum(id,x,t);
			break;
		default :
			std::cout << "undefined ensemble type" << std::endl;
			return 0;
			break;
	}
}

prec Dinamic::to_odeint::p_force_sum(prec sum, int i,const std::vector<prec> &x,const prec t)
{
	return sum/I_t_odeint[i]+F_t_odeint[i]*sin(W_t_odeint[i]*t-x[2*i])/I_t_odeint[i]-(G_t_odeint[i]/I_t_odeint[i])*x[2*i+1];
}

prec Dinamic::to_odeint::h_force_sum(prec sum, int i,const std::vector<prec> &x,const prec t)
{
	return sum/I_t_odeint[i]+F_t_odeint[i]*(sin(W_t_odeint[i]*t)-x[2*i])/I_t_odeint[i]-(G_t_odeint[i]/I_t_odeint[i])*x[2*i+1];
}

prec Dinamic::to_odeint::t_force_sum(prec sum, int i,const std::vector<prec> &x,const prec t)
{
	if(F_t_odeint[i]<0.0000001 && F_t_odeint[i]>-0.0000001)
	{
		return sum/I_t_odeint[i]-(G_t_odeint[i]/I_t_odeint[i])*x[2*i+1];	
	}
	else
	{
		int n=static_cast<int>(t*W_t_odeint[i]/(2*M_PI));
		prec fixed_t=t-n*2.0*M_PI/W_t_odeint[i];

		if(fixed_t<1.0*M_PI/W_t_odeint[i])
		{
			return sum/I_t_odeint[i]+F_t_odeint[i]*(-1+2.0*fixed_t/(1.0*M_PI/W_t_odeint[i]) -x[2*i])/I_t_odeint[i]-(G_t_odeint[i]/I_t_odeint[i])*x[2*i+1];
		}
		else
		{
			return sum/I_t_odeint[i]+F_t_odeint[i]*(1-2.0*(fixed_t-1.0*M_PI/W_t_odeint[i])/(1.0*M_PI/W_t_odeint[i]) -x[2*i])/I_t_odeint[i]-(G_t_odeint[i]/I_t_odeint[i])*x[2*i+1];
		}
	}
}

prec Dinamic::to_odeint::s_force_sum(prec sum, int i,const std::vector<prec> &x,const prec t)
{
	if(F_t_odeint[i]<0.0000001 && F_t_odeint[i]>-0.0000001)
	{
		return sum/I_t_odeint[i]-(G_t_odeint[i]/I_t_odeint[i])*x[2*i+1];
	}
	else
	{
		int n=static_cast<int>(t*W_t_odeint[i]/(2*M_PI));
		prec fixed_t=t-n*2.0*M_PI/W_t_odeint[i];

		return sum/I_t_odeint[i]+F_t_odeint[i]*(-1+2.0*fixed_t/(2.0*M_PI/W_t_odeint[i]) -x[2*i])/I_t_odeint[i]-(G_t_odeint[i]/I_t_odeint[i])*x[2*i+1];
	}
}

prec Dinamic::to_odeint::c_a_force_sum(prec sum, int i,const std::vector<prec> &x,const prec t)
{
//return sum/I_t_odeint[i]+F_t_odeint[i]*(sin(W_t_odeint[i]*t)+(1.0/5.0)*sin(2*W_t_odeint[i]*t)-(1.0/5.0)*sin(3*W_t_odeint[i]*t)-x[2*i])/I_t_odeint[i]-(G_t_odeint[i]/I_t_odeint[i])*x[2*i+1];
	//F.assign(0,25);
	return sum/I_t_odeint[i]+F_t_odeint[i]*(sin(W_t_odeint[i]*t)+(0.473916/2.0)*sin(2*W_t_odeint[i]*t)+(0.221962/3.0)*sin(3*W_t_odeint[i]*t)+(0.102042/4.0)*sin(4*W_t_odeint[i]*t)-x[2*i])/I_t_odeint[i]-(G_t_odeint[i]/I_t_odeint[i])*x[2*i+1];
}

prec Dinamic::to_odeint::solid_force_sum(int i,const std::vector<prec> &x,const prec t)
{
	return -K_t_odeint*sin(x[2*i])/I_t_odeint[i]+F_t_odeint[i]*sin(W_t_odeint[i]*t-x[2*i])/I_t_odeint[i]-(G_t_odeint[i]/I_t_odeint[i])*x[2*i+1];
}


prec Dinamic::to_odeint::interaction_sum_all(int id, const std::vector<prec> &x)
{
	prec sum=0;
	for (int i = 0; i < N_t_odeint; ++i)
	{
		sum=sum+A_t_odeint(id,i)*sin(x[2*i]-x[2*id])/N_t_odeint;
	}
	return sum*K_t_odeint;
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

void Dinamic::init_I_type_rand(bool c_bool, bool c_f_belong)
{
	I.type_rand(c_bool, c_f_belong);
	I.generate();
}

void Dinamic::init_F_type_rand(bool c_bool, bool c_f_belong)
{
	F.type_rand(c_bool, c_f_belong);
	F.generate();
}

void Dinamic::init_G_type_rand(bool c_bool, bool c_f_belong)
{
	G.type_rand(c_bool, c_f_belong);
	G.generate();
}

void Dinamic::init_W_type_rand(bool c_bool, bool c_f_belong)
{
	W.type_rand(c_bool, c_f_belong);
	W.generate();
}

std::string Dinamic::get_type()
{
	return A.get_type();
}

void Dinamic::new_va_F(int id,prec val)
{
	F.assign(id,val);
}
