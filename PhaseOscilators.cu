#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>
#include <limits>

#include <sys/stat.h>
#include <boost/numeric/odeint.hpp>
#include <boost/random.hpp>

typedef double prec;

#include "Ensemble_param.h"
#include "Ensemble_connections.h"
#include "Ensemble_Dinamic.h"

class Evolve
{
	int N;
	Ensemble_param  x_i;
	Ensemble_param  v_i;
	std::vector<std::vector<prec>> q;
	std::vector<prec> t;

	std::vector<prec> translated_init;
	size_t steps=0;
	

	void translate_init()
	{
		for (int i = 0; i < N; ++i)
		{
			translated_init[2*i]=x_i[i];
			translated_init[2*i+1]=v_i[i];
		}
	}

	struct push_back_state_and_time
	{
	    std::vector< std::vector<prec> >& m_states;
	    std::vector< prec >& m_times;

	    prec endtime;
	
	    push_back_state_and_time( std::vector< std::vector<prec> > &states , std::vector< prec > &times, prec in_endtime ) : m_states( states ) , m_times( times ), endtime(in_endtime) { }
	
	    void operator()( const std::vector<prec> &x , prec t_boost )
	    {
	    	if(t_boost>endtime-500)
	    	{
	        	m_states.push_back( x );
	        	m_times.push_back( t_boost );
	    	}
	    }
	};


public:
	Evolve(int in_N, boost::mt19937 &in_rng) : N(in_N), x_i(N,"X_i",0,0,in_rng), v_i(N,"V_i",0,0,in_rng), translated_init(2*N)
	{
		translate_init();
	}

	void run(Dinamic &P)
	{
		boost::numeric::odeint::runge_kutta4 < std::vector< double > > stepper;
		steps=integrate_adaptive(stepper, P, translated_init , 1000.0 , 1500.0 , 0.01,push_back_state_and_time( q , t,500));
	}

	void print_init()
	{
		x_i.print();
		v_i.print();
	}

	void load_init()
	{
		x_i.load();
		v_i.load();
		translate_init();
	}

	void reset_init()
	{
		for (int i = 0; i < N; ++i)
		{
			x_i.assign(i,x(i,q.size()-1));
			v_i.assign(i,v(i,q.size()-1));
		}
	}

	prec x(int i,int t_i)
	{
		return q[t_i][2*i];
	}

	prec v(int i,int t_i)
	{
		return q[t_i][2*i+1];
	}

	void print_x()
	{
		mkdir("params", 0777);
	
		std::ofstream txtOut;
		txtOut.precision(std::numeric_limits< prec >::max_digits10);
		txtOut.open("params/X.txt");
		for (int j = 0; j < steps; ++j)
		{
			txtOut << t[j] << " ";
			for (int i = 0; i < N; ++i)
			{
				txtOut << x(i,j) << " ";
			}
			txtOut << std::endl;
		}

		txtOut.close();
	}

	void print_v()
	{
		mkdir("params", 0777);
	
		std::ofstream txtOut;
		txtOut.precision(std::numeric_limits< prec >::max_digits10);
		txtOut.open("params/V.txt");
		for (int j = 0; j < steps; ++j)
		{
			txtOut << t[j] << " ";
			for (int i = 0; i < N; ++i)
			{
				txtOut << v(i,j) << " ";
			}
			txtOut << std::endl;
		}

		txtOut.close();
	}

	void print()
	{
		print_x();
		print_v();
	}



};



int main()
{
 	boost::mt19937 rng(static_cast<unsigned int>(std::time(0)));

 	int N=10;

	Dinamic P(N ,1 ,rng ,1 ,0 ,1000 ,0 ,2.5 ,0 ,1 ,0,"chain"); //N, K, rng, I, sigm_I, F, sigm_F, G, sigm_G, W, sigm_W
	P.print_params();
	P.print_params_to_console();
	Evolve e(N, rng);
	e.load_init();
	e.run(P);
	e.reset_init();
	e.print_init();
	e.print();

	return 0;
}