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


	boost::numeric::odeint::runge_kutta4 <std::vector< prec > , prec , std::vector< prec > , prec , boost::numeric::odeint::openmp_range_algebra> stepper;
	push_back_state_and_time push(q , t, t_start, t_end, t_save);
	size_t steps=integrate_adaptive(stepper, P(), translated_init , t_start , t_end, 0.01 , push);
	std::cout << std::endl;
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