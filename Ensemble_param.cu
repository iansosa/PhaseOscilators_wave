#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>

#include <sys/stat.h>
#include <boost/numeric/odeint.hpp>
#include <boost/random.hpp>

typedef double prec;

#include "Ensemble_param.h"

void Ensemble_param::equal(std::vector<prec> &x)
{
	container=x;
}

void Ensemble_param::assign(int i, prec value)
{
	container[i]=value;
}

void Ensemble_param::dist(prec in_mid, prec in_sigma)
{
	mid=in_mid;
	sigma=in_sigma;
}

int Ensemble_param::size()
{
	return N;
}

void Ensemble_param::resize(int in_N)
{
	N=in_N;
	container.resize(in_N);
	generate(false);
}

void Ensemble_param::generate(bool first)
{
	boost::normal_distribution<> gauss(mid,sigma);
    boost::variate_generator< boost::mt19937&, boost::normal_distribution<> > gen(rng,gauss);
    for (int i = 0; i < container.size(); ++i)
    {
    	container[i]=gen();
    }

    if(first==false)
    {
    	container[0]=mid;
    }

    if(type=="ofirst")
    {
    	container[0]=mid;
    	for (int i = 1; i < container.size(); ++i)
    	{
    		container[i]=0;
    	}
    }
}

void Ensemble_param::print()
{
	mkdir("params", 0777);

	std::ofstream txtOut;
	txtOut.open("params/"+name+".txt");
	txtOut.precision(std::numeric_limits< prec >::max_digits10);
	txtOut << N << " ";
	txtOut << mid << " ";
	txtOut << sigma << " ";
	txtOut << type << " ";
	for (int i = 0; i < container.size(); ++i)
	{
		txtOut << container[i] << " ";
	}
	txtOut.close();
}

void Ensemble_param::load()
{
	std::ifstream txtIn_check;
	int N_check;
	txtIn_check.open("params.txt");
	txtIn_check >> N_check;
	txtIn_check.close();	
	std::ifstream txtIn;
	txtIn.open("params/"+name+".txt");
	txtIn >> N;
	txtIn >> mid;
	txtIn >> sigma;
	txtIn >> type;
	container.resize(N);
	if(N_check!=N)
	{
		std::cout << "INCONSISTENT PARAM SIZE FOR " << name << ". Expects N=" << N << " but check throws N=" << N_check <<std::endl;
	}
	for (int i = 0; i < container.size(); ++i)
	{
		txtIn >> container[i];
	}
	txtIn.close();
}

void Ensemble_param::type_rand(bool is_rand)
{
	if(is_rand==true)
	{
		type="gauss";
	}
	else
	{
		type="ofirst";
	}
}

void Ensemble_param::print_to_console()
{
	std::cout << name <<": N="<< N << " mid="<< mid << " sigma=" << sigma << " " << type<<std::endl;
	std::cout << "   ";
	for (int i = 0; i < container.size(); ++i)
	{
		std::cout << container[i] << " ";
	}
	std::cout << std::endl;
}

prec& Ensemble_param::operator[](int i)
{
	return container[i];
}