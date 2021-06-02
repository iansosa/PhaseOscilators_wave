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
#include "Ensemble_evolve.h"

void Evolve::translate_init()
{
	for (int i = 0; i < N; ++i)
	{
		translated_init[2*i]=x_i[i];
		translated_init[2*i+1]=v_i[i];
		if(i==0)
		{
			int entero=translated_init[2*i]/(2.0*M_PI);
			translated_init[2*i]=translated_init[2*i]-2*M_PI*entero;
		}
	}
}

void Evolve::run(Dinamic &P, prec t_start, prec t_end, prec t_save)
{
	if(t_save<0)
	{
		t_save=t_end;
	}
	std::cout << std::endl;
	std::cout << "N=" << N << std::endl;
	std::cout << "Running..." << std::endl;

	prec dt=0.01;
	boost::numeric::odeint::runge_kutta4 <std::vector< prec > , prec , std::vector< prec > , prec , boost::numeric::odeint::openmp_range_algebra> stepper;
	push_back_state_and_time push(q , t, t_start, t_end, t_save,dt);
	
	size_t steps=integrate_adaptive(stepper, P(), translated_init , t_start , t_end, dt , push);
	std::cout << std::endl;
}

int Evolve::find_next_maxima(int t_start,int k)
{
	while(v(k,t_start+1)<=v(k,t_start) && t_start < t.size()-2)
	{
		t_start++;
	}
	if(v(k,t_start+1)>v(k,t_start))
	{
		for (int i = t_start; i < t.size()-1; ++i)
		{
			if(v(k,i+1)<=v(k,i))
			{
				return i;
			}
		}
	}
	return t.size()-1;
}

int Evolve::find_next_minima(int t_start,int k)
{
	while(v(k,t_start+1)>=v(k,t_start) && t_start < t.size()-2)
	{
		t_start++;
	}
	if(v(k,t_start+1)<v(k,t_start))
	{
		for (int i = t_start; i < t.size()-1; ++i)
		{
			if(v(k,i+1)>=v(k,i))
			{
				return i;
			}
		}
	}
	return t.size()-1;
}

bool Evolve::did_converge(int i)
{
	if(convergence(i)>conv_crit)
	{
		return true;
	}
	else
	{
		return false;
	}
}

prec Evolve::convergence(int k)
{
	int first_maxima=find_next_maxima(0,k);
	if(first_maxima<t.size()-1)
	{
		int second_maxima=find_next_maxima(first_maxima+1,k);
		int first_minima=find_next_minima(first_maxima+1,k);
		if(second_maxima<t.size()-1 && first_minima<t.size()-1 && second_maxima-first_minima>50)
		{
			prec first_v=v(k,first_maxima);
			prec second_v=v(k,first_minima);
			if(first_v>second_v)
			{
				return 100-100*fabs((v(k,second_maxima)-v(k,first_maxima))/(first_v-second_v));
			}
		}
	}
	return 0;
}

prec Evolve::period(int k)
{
	if(convergence(k)>conv_crit)
	{
		int first_maxima=find_next_maxima(0,k);
		int second_maxima=find_next_maxima(first_maxima+1,k);
		return t[second_maxima]-t[first_maxima];
	}
	return 0;
}

prec Evolve::frec(int k)
{
	prec ret=period(k);
	if(ret>0)
	{
		return 2*M_PI/ret;
	}
	return 0;
}

prec Evolve::drift(int k)
{
	//std::cout << k << " entered w" << std::endl;
	if(convergence(k)>conv_crit)
	{
		//std::cout << "drift: " << k << " passed conv test" << std::endl;
		int first_maxima=find_next_maxima(1,k);
		//std::cout <<"drift: K:" << k << " first_maxima: " << first_maxima;
		if(first_maxima<t.size()-1)
		{
			int last_maxima=find_next_maxima(first_maxima+1,k);
			
			int moving=last_maxima;
			
			while(moving<t.size()-1)
			{
				last_maxima=moving;
				moving=find_next_maxima(moving+1,k);
			}
			//std::cout << "last_maxima: " << last_maxima << std::endl;
			if(first_maxima==last_maxima)
			{
				return 0;
			}
			return (x(k,last_maxima)-x(k,first_maxima))/(t[last_maxima]-t[first_maxima]);
		}
	}
	return 0;
}

prec Evolve::Amp(int k,Dinamic &P)
{
	if(did_converge(k)==true)
	{
		//std::cout << k << " did converge" << std::endl;
		int tstart=find_next_maxima(0,k);
		//std::cout << "amp: first maxima " << tstart << " v " << v(k,tstart);
		int tend=find_next_maxima(tstart+1,k);
		//std::cout << " second maxima " << tend << " v" << v(k,tend) << std::endl;
		if(tstart==tend)
		{
			return 0;
		}
		prec maxima=-10000;
		prec minima=100000;
		prec current;
		prec w;
		if(P.get_type()=="h_chain")
		{
			w=0;
		}
		else
		{
			w=drift(k);
		}
		for (int i = tstart; i < tend; ++i)
		{
			
			current=x(k,i)-w*t[i];
			if(maxima<current)
			{
				maxima=current;
			}
			if(minima>current)
			{
				minima=current;
			}
		}
		return (maxima-minima)/2.0;
	}
	return 0;
}

prec Evolve::Diff(int k)
{
	if(did_converge(k)==true && k<N-1 && k>0)
	{
		int tstart=find_next_maxima(0,k);
		int tend=find_next_maxima(tstart+1,k);
		if(tstart==tend)
		{
			return 0;
		}
		prec maxima=-10000;
		prec current;
		for (int i = tstart; i < tend; ++i)
		{
			current=fabs(x(k,i)-x(k+1,i));
			if(maxima<current)
			{
				maxima=current;
			}
		}
		return maxima;
	}
	return 0;
}

void Evolve::print_init()
{
	x_i.print();
	v_i.print();
}

void Evolve::load_init()
{
	x_i.load();
	v_i.load();
	translate_init();
}

void Evolve::reset_init()
{
	for (int i = 0; i < N; ++i)
	{
		x_i.assign(i,x(i,q.size()-1));
		v_i.assign(i,v(i,q.size()-1));
	}
	translate_init();
}

prec Evolve::x(int i,int t_i)
{
	return q[t_i][2*i];
}

prec Evolve::v(int i,int t_i)
{
	return q[t_i][2*i+1];
}

void Evolve::print_x()
{
	mkdir("out", 0777);

	std::ofstream txtOut;
	txtOut.precision(std::numeric_limits< prec >::max_digits10);
	txtOut.open("out/X.txt");
	for (int j = 0; j < t.size(); ++j)
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

void Evolve::print_v()
{
	mkdir("out", 0777);

	std::ofstream txtOut;
	txtOut.precision(std::numeric_limits< prec >::max_digits10);
	txtOut.open("out/V.txt");
	for (int j = 0; j < t.size(); ++j)
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

void Evolve::print()
{
	std::cout << "Printing X..." << std::endl;
	print_x();
	std::cout << "Printing V..." << std::endl;
	print_v();
}

void Evolve::clean()
{
	q.clear();
	t.clear();

}

int Evolve::size()
{
	return N;
}