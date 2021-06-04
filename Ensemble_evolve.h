#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#ifndef ENSEMBLE_EVOLVE_H
#define ENSEMBLE_EVOLVE_H

class Evolve
{
	int N;
	Ensemble_param  x_i;
	Ensemble_param  v_i;
	std::vector<std::vector<prec>> q;
	std::vector<prec> t;
	
	prec conv_crit=99.7;

	std::vector<prec> translated_init;

	void translate_init();

	struct push_back_state_and_time
	{
	    std::vector< std::vector<prec> >& m_states;
	    std::vector< prec >& m_times;

	    prec t_start;
	    prec t_end;
	    prec t_save;
	    prec dt;
	    int iter=0;
	
	    push_back_state_and_time( std::vector< std::vector<prec> > &states , std::vector< prec > &times, prec in_t_start, prec in_t_end , prec in_t_save, prec dt_in) : m_states( states ) , m_times( times ), t_start(in_t_start), t_end(in_t_end), t_save(in_t_save), dt(dt_in) { }
	
	    void operator()( const std::vector<prec> &x , prec t_boost )
	    {
	    	if(iter%10000==0)
	    	{
	    		std::cout << "time: " << iter*dt+t_start << "/" << t_end <<'\r' << std::flush ;
	    	}
	    	iter++;
	    	if(t_boost>=t_end-t_save)
	    	{
	        	m_states.push_back( x );
	        	m_times.push_back( t_boost );
	    	}
	    }
	};

public:
	std::vector<bool> conv;

	Evolve(int in_N, boost::mt19937 &in_rng) : N(in_N), x_i(in_N,"X_i",0,0,in_rng), v_i(in_N,"V_i",0,0,in_rng), translated_init(2*in_N), conv(in_N)
	{
		translate_init();
	}

	void run(Dinamic &P, prec t_start, prec t_end, prec t_save=-1);

	int find_next_maxima(int t_start,int k=1);
	int find_next_minima(int t_start,int k=1);

	void calc_convergence();

	prec convergence(int i);

	prec period(int i=1);
	prec frec(int i=1);
	prec drift(int i=1);

	void print_init();

	void load_init();

	void reset_init();

	prec x(int i,int t_i);

	prec v(int i,int t_i);

	void print_x();

	void print_v();

	void print();

	void clean();

	int size();

	prec Amp(int i,Dinamic &P);

	prec Diff(int i);
};

#endif