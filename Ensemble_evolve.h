#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#ifndef ENSEMBLE_EVOLVE_H
#define ENSEMBLE_EVOLVE_H

class Evolve
{
	int N; //number of oscilators
	Ensemble_param  x_i; //untranslated initial conditions
	Ensemble_param  v_i; //untranslated initial conditions
	std::vector<std::vector<prec>> q; //position and velocity vector
	std::vector<prec> t; //time vector
	prec conv_crit=99.9; //convergence criterium
	int max_stats_prom=30;
	int max_RMS_prom=20;
	std::vector<prec> translated_init; // vector of initial conditions of space and velocity
	std::vector<bool> conv_method;

	void translate_init(); //copies x_i and v_i consistently to translated_init
	int find_next_maxima(int t_start,int k=1); //finds next maxima for the time evolution from time t_start and oscilator k
	int find_next_minima(int t_start,int k=1); //finds next minima for the time evolution from time t_start and oscilator k
	prec convergence(int i); //determines convergence number of oscilator i
	void print_x(bool coherent); //prints positions to file
	void print_v(bool coherent); //prints velocities to file

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

	Evolve(int in_N, boost::mt19937 &in_rng) : N(in_N), x_i(in_N,"X_i",0,0,in_rng), v_i(in_N,"V_i",0,0,in_rng), translated_init(2*in_N), conv(in_N)
	{
		translate_init();

		conv_method.resize(N,true);
	}

	std::vector<bool> conv; //vector of convergence states

	void run(Dinamic &P, prec t_start, prec t_end, prec t_save=-1,prec dt=0.01); //runs simulation with start time, end time, and save time
	void calc_convergence(); //determines convergence state vector
	void set_conv_method(int k);
	prec period(int i=1); //determines period of oscilator i
	prec frec(int i=1); //determines frecuency of oscilator i
	prec drift(int i=1); //determines drift of oscilator i
	void print_init(); //prints initial conditions x and v
	void load_init(); //loads initial condicions to x and v from file and translates it
	void set_init(std::vector<prec> x_in,std::vector<prec> v_in); //sets initial conditions from vectors
	void reset_init(); //resets initial conditions to last time step in the evolution
	prec x(int i,int t_i); //returns the position of oscilator i at timestep t_i
	prec v(int i,int t_i); //returns the velocity of oscilator i at timestep t_i
	prec ts(int t_i); //returns the time value at timestep t_i
	int t_size(); //returns size of time vector
	void print(bool coherent=false); //prints positions and velocities to files
	void clean(); //clears state vector and time vector
	int size(); //return number of oscilators
	prec Amp(int i,Dinamic &P); //returns amplitude of oscilator i in the dinamic P
	prec RMS(int i,Dinamic &P); //returns RMS of oscilator i in the dinamic P
	prec Diff_max(int i); //returns maximun difference between oscilator k and k+1 from the start of saved time
	prec Diff_med(int i); //returns average difference between oscilator k and k+1 from the start of saved time
	prec Diff_min(int i); //returns minimum difference between oscilator k and k+1 from the start of saved time
};

#endif