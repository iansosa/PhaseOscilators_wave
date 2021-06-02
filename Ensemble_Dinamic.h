#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <ctime>

#include <boost/numeric/odeint/external/openmp/openmp.hpp>

#ifndef ENSEMBLE_DINAMIC_H
#define ENSEMBLE_DINAMIC_H

class Dinamic
{
	int N;
    prec K;

    Ensemble_param I;
    Ensemble_param F;
    Ensemble_connections A;
    Ensemble_param G;
    Ensemble_param W;

    boost::mt19937 &rng;

    public:

    	class to_odeint
    	{
			int &N_t_odeint;
    		prec &K_t_odeint;
	
    		Ensemble_param &I_t_odeint;
    		Ensemble_param &F_t_odeint;
    		Ensemble_connections &A_t_odeint;
    		Ensemble_param &G_t_odeint;
   	 		Ensemble_param &W_t_odeint;
   	 	public:
   	 		to_odeint(int &N_in_t_odeint, prec &K_in_t_odeint,Ensemble_param &in_I_t_odeint,Ensemble_param &in_F_t_odeint,Ensemble_connections &in_A_t_odeint,Ensemble_param &in_G_t_odeint,Ensemble_param &in_W_t_odeint) : N_t_odeint(N_in_t_odeint), K_t_odeint(K_in_t_odeint),I_t_odeint(in_I_t_odeint),F_t_odeint(in_F_t_odeint),A_t_odeint(in_A_t_odeint),G_t_odeint(in_G_t_odeint),W_t_odeint(in_W_t_odeint)
   	 		{}

			prec interaction_sum_chain(int id,const std::vector<prec> &x);
			prec h_interaction_sum_chain(int id,const std::vector<prec> &x);
			prec interaction_sum(int id, const std::vector<prec> &x);
			prec force_sum(prec sum,int id, const std::vector<prec> &x,const prec t);
			prec h_force_sum(prec sum,int id, const std::vector<prec> &x,const prec t);
			prec p_force_sum(prec sum,int id, const std::vector<prec> &x,const prec t);
			prec interaction_sum_all(int id, const std::vector<prec> &x);

    		void operator() (const std::vector<prec> &x ,std::vector<prec> &dxdt ,const prec t)
    		{
    			#pragma omp parallel for schedule(runtime)
    	    	for (int i = 0; i < N_t_odeint; ++i)
    		    {
    	 	   		prec sum;
    		    	sum=interaction_sum(i,x);
        			dxdt[2*i]=x[2*i+1];
        			dxdt[2*i+1]=force_sum(sum,i,x,t);
					//dxdt[2*i+1]=sum/I_t_odeint[i]+F_t_odeint[i]*sin(W_t_odeint[i]*t-x[2*i])/I_t_odeint[i]-(G_t_odeint[i]/I_t_odeint[i])*x[2*i+1];
        		}
    		}

    	};

	
    	Dinamic(boost::mt19937 &in_rng) : rng(in_rng), A("A", in_rng),  I("I",in_rng), F("F",in_rng), G("G",in_rng), W("W",in_rng)
    	{
			std::ifstream txtIn_check;
			txtIn_check.open("params.txt");
			txtIn_check >> N;
			txtIn_check >> K;
			txtIn_check.close();
    	}

    	Dinamic(int in_N , prec in_K, boost::mt19937 &in_rng, prec mid_I, prec sigma_I, prec mid_F, prec sigma_F, prec mid_G, prec sigma_G, prec mid_W, prec sigma_W,std::string A_type="global") : N(in_N), K(in_K), rng(in_rng), A(in_N,"A",rng,A_type), I(in_N,"I",mid_I,sigma_I,in_rng), F(in_N,"F",mid_F,sigma_F,in_rng,false), G(in_N,"G",mid_G,sigma_G,in_rng), W(in_N,"W",mid_W,sigma_W,in_rng,false)
    	{ }

    	to_odeint operator() ()
    	{
    		to_odeint a(N,K,I,F,A,G,W);
    		return a;
    	}

    	void print_params();

    	void print_params_to_console();

    	void generate();

    	std::string get_type();

};

#endif