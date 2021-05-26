#include <string>
#include <vector>
#include <iostream>
#include <fstream>

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

    prec interaction_sum_all(int id, const std::vector<prec> &x);

    prec interaction_sum_chain(int id,const std::vector<prec> &x);

	public:

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

    	void operator() (const std::vector<prec> &x ,std::vector<prec> &dxdt ,const prec t)
    	{
    	    for (int i = 0; i < N; ++i)
    	    {
    	    	prec sum=interaction_sum(i,x);

        		dxdt[2*i]=x[2*i+1];
				dxdt[2*i+1]=sum/I[i]+F[i]*sin(W[i]*t-x[2*i])/I[i]-(G[i]/I[i])*x[2*i+1];      	
        	}
    	}

    	prec interaction_sum(int id, const std::vector<prec> &x);

    	void print_params();

    	void print_params_to_console();

    	void generate();

};

#endif