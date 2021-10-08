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
#include "Ensemble_evolve.h"

void Evolve::translate_init()
{
	//std::cout << "entered translate_init" << std::endl;
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

int Evolve::find_next_maxima(int t_start,int k)
{
	//std::cout << k <<" entered find_next_maxima t_start " << t_start << std::endl;

	int save_start=t_start;
	if(conv_method[k]==true)
	{
		while(v(k,t_start+1)<=v(k,t_start) && t_start < t.size()-3)
		{
			t_start=t_start+1;
		}
		if(v(k,t_start+1)>v(k,t_start))
		{
			for (int i = t_start; i < t.size()-3; ++i)
			{
				if(v(k,i+1)<=v(k,i))
				{
					return i;
				}
			}
		}
		return save_start;
	}
	else
	{
		while(x(k,t_start+1)<=x(k,t_start) && t_start < t.size()-3)
		{
			t_start=t_start+1;
		}
		if(x(k,t_start+1)>x(k,t_start))
		{
			for (int i = t_start; i < t.size()-3; ++i)
			{
				if(x(k,i+1)<=x(k,i))
				{
					//std::cout << k << " " << i << std::endl;
					return i;
				}
			}
		}
		return save_start;
	}
}

int Evolve::find_next_minima(int t_start,int k)
{
	int save_start=t_start;
	if(conv_method[k]==true)
	{
		while(v(k,t_start+1)>=v(k,t_start) && t_start < t.size()-3)
		{
			t_start=t_start+1;
		}
		if(v(k,t_start+1)<v(k,t_start))
		{
			for (int i = t_start; i < t.size()-3; ++i)
			{
				if(v(k,i+1)>=v(k,i))
				{
					return i;
				}
			}
		}
		return save_start;
	}
	else
	{
		while(x(k,t_start+1)>=x(k,t_start) && t_start < t.size()-3)
		{
			t_start=t_start+1;
		}
		if(x(k,t_start+1)<x(k,t_start))
		{
			for (int i = t_start; i < t.size()-3; ++i)
			{
				if(x(k,i+1)>=x(k,i))
				{
					return i;
				}
			}
		}
		return save_start;
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

void Evolve::print_x(bool coherent)
{
	int t_end=t.size();
	int t_start=0;
	if(coherent==true)
	{
		t_start=find_next_maxima(t_start,10);
		t_end=1;
		for (int i = 0; i < max_stats_prom; ++i)
		{
			t_end=find_next_maxima(t_end,10);
		}
	}
	mkdir("out", 0777);

	std::ofstream txtOut;
	txtOut.precision(std::numeric_limits< prec >::max_digits10);
	txtOut.open("out/X.txt");
	for (int j = t_start; j < t_end; ++j)
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

void Evolve::print_v(bool coherent)
{
	int t_end=t.size();
	int t_start=0;
	if(coherent==true)
	{
		t_start=find_next_maxima(t_start,10);
		t_end=1;
		for (int i = 0; i < max_stats_prom; ++i)
		{
			t_end=find_next_maxima(t_end,10);
		}
	}
	mkdir("out", 0777);

	std::ofstream txtOut;
	txtOut.precision(std::numeric_limits< prec >::max_digits10);
	txtOut.open("out/V.txt");
	for (int j = t_start; j < t_end; ++j)
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

void Evolve::run(Dinamic &P, prec t_start, prec t_end, prec t_save,prec dt)
{
	if(t_save<0)
	{
		t_save=t_end;
	}
	std::cout << std::endl;
	std::cout << "N=" << N << std::endl;
	std::cout << "Running..." << std::endl;

	boost::numeric::odeint::runge_kutta4 <std::vector< prec > , prec , std::vector< prec > , prec , boost::numeric::odeint::openmp_range_algebra> stepper;
	push_back_state_and_time push(q , t, t_start, t_end, t_save,dt);
	
	size_t steps=integrate_adaptive(stepper, P(), translated_init , t_start , t_end, dt , push);
	std::cout << std::endl;
}

void Evolve::calc_convergence()
{
	//std::cout << "entered calc_convergence" << std::endl;
	conv[0]=false;
	for (int i = 1; i < N; ++i)
	{
		if(convergence(i)>conv_crit)
		{
			conv[i]=true;
		}
		else
		{
			conv[i]=false;
		}
		if(conv_method[i]==false)
		{
			conv[i]=true;
		}
	}
}

void Evolve::set_conv_method(int k)
{
	for (int i = 0; i < k; ++i)
	{
		conv_method[i]=false;
	}
}

prec Evolve::period(int k)
{
	if(conv[k]==true)
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
	if(conv[k]==true)
	{
		//std::cout << "drift: " << k << " passed conv test" << std::endl;
		int first_maxima=find_next_maxima(1,k);
		//std::cout <<"drift: K:" << k << " first_maxima: " << first_maxima << std::endl;
		if(first_maxima<t.size()-1)
		{
			int last_maxima=find_next_maxima(first_maxima,k);
			
			int moving=last_maxima;
			
			for (int i = 0; i < max_stats_prom; ++i)
			{
				last_maxima=find_next_maxima(last_maxima,k);
			}
			/*while(moving<t.size()-10)
			{
				//std::cout << t.size() << std::endl;
				last_maxima=moving;
				moving=find_next_maxima(moving,k);
			}*/
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

void Evolve::set_init(std::vector<prec> x_in,std::vector<prec> v_in)
{
	if(x_i.size() != x_in.size())
	{
		std::cout << "Size difference in initial vector encountered" << std::endl;
		return;
	}

	for (int i = 0; i < x_i.size(); ++i)
	{
		x_i.assign(i,x_in[i]);
		v_i.assign(i,v_in[i]);
	}
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

prec Evolve::ts(int t_i)
{
	return t[t_i];
}

int Evolve::t_size()
{
	return t.size();
}

void Evolve::print(bool coherent)
{
	std::cout << "Printing X..." << std::endl;
	print_x(coherent);
	std::cout << "Printing V..." << std::endl;
	print_v(coherent);
}

void Evolve::clean()
{
	std::vector<std::vector<prec>> ().swap(q);
	std::vector<prec> ().swap(t);
	//q.clear();
	//t.clear();
	translate_init();
}

int Evolve::size()
{
	return N;
}

prec Evolve::Amp(int k,Dinamic &P)
{
	if(conv[k]==true)
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
		if(P.get_type()=="h_chain" || P.get_type()=="t_chain" || P.get_type()=="s_chain" || P.get_type()=="c_a_chain")
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

prec Evolve::RMS(int k,Dinamic &P)
{
	if(conv[k]==true)
	{
		//std::cout << k << " did converge" << std::endl;
		int tstart=find_next_maxima(0,k);
		//std::cout << "amp: first maxima " << tstart << " v " << v(k,tstart);
		int tend=find_next_maxima(tstart+1,k);

		for (int i = 0; i < max_RMS_prom; ++i)
		{
			tend=find_next_maxima(tend,k);
		}
		//std::cout << " second maxima " << tend << " v" << v(k,tend) << std::endl;
		if(tstart==tend)
		{
			return 0;
		}
		
		prec w;
		if(P.get_type()=="h_chain" || P.get_type()=="t_chain" || P.get_type()=="s_chain" || P.get_type()=="c_a_chain")
		{
			w=0;
		}
		else
		{
			w=drift(k);
		}

		prec base=0;
		prec current;
		for (int i = tstart; i <= tend; ++i)
		{
			current=x(k,i)-w*t[i];
			base=base+current;
		}
		//std::cout.precision(std::numeric_limits< prec >::max_digits10);
		base=base/(tend+1-tstart);
		//std::cout << k << " tstart=" << tstart << " tstart_value=" << x(k,tstart)-w*t[tstart] << " xvalue=" << x(k,tstart) << " wvalue=" << w << " tvalue=" << t[tstart] << " tend=" << tend << " tend_value=" << x(k,tend)-w*t[tend] <<" base=" << base << std::endl;
		
		prec sum=0;
		for (int i = tstart; i <= tend; ++i)
		{
			sum=sum+pow(x(k,i)-w*t[i]-base,2)*0.01;
		}
		sum=sum/(t[tend]-t[tstart]);
		sum=sqrt(sum);
		return sum;
	}
	return 0;
}

prec Evolve::Diff_max(int k)
{
	if(conv[k]==true && k<N-1 && k>0)
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

prec Evolve::Diff_med(int k)
{
	if(conv[k]==true && k<N-1 && k>0)
	{
		int tstart=find_next_maxima(0,k);
		int tend=find_next_maxima(tstart+1,k);
		if(tstart==tend)
		{
			return 0;
		}
		prec current=0;
		int count=0;
		for (int i = tstart; i <= tend; ++i)
		{
			count++;
			current=current+(x(k,i)-x(k+1,i));
		}
		return current/count;
	}
	return 0;
}

prec Evolve::Diff_min(int k)
{
	if(conv[k]==true && k<N-1 && k>0)
	{
		int tstart=find_next_maxima(0,k);
		int tend=find_next_maxima(tstart+1,k);
		if(tstart==tend)
		{
			return 0;
		}
		prec minima=10000;
		prec current;
		for (int i = tstart; i < tend; ++i)
		{
			current=fabs(x(k,i)-x(k+1,i));
			if(minima>current)
			{
				minima=current;
			}
		}
		return minima;
	}
	return 0;
}