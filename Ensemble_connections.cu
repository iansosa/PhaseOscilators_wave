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

#include "Ensemble_connections.h"


int Ensemble_connections::_2dto1d(int a, int b)
{
	return a+b*N;
}

void Ensemble_connections::_1dto2d(int idx,int &a, int &b)
{
	a=idx%N;
	b=idx/N;
}

int Ensemble_connections::size()
{
	return N;
}

void Ensemble_connections::resize(int i)
{
	N=i;
	container.resize(N*N);
	generate(type);
}

void Ensemble_connections::generate(std::string s)
{
	if(s=="global")
	{
		std::fill(container.begin(), container.end(), 1);
	}
	if(s=="empty")
	{
		std::fill(container.begin(), container.end(), 0);
	}
	if(s=="chain")
	{
		std::fill(container.begin(), container.end(), 0);
		for (int i = 0; i < N-1; ++i)
		{
			connect(i,i+1);
		}
	}
	if(s=="h_chain")
	{
		std::fill(container.begin(), container.end(), 0);
		for (int i = 0; i < N-1; ++i)
		{
			connect(i,i+1);
		}
	}
	
}

void Ensemble_connections::connect(int i, int j)
{
	container[_2dto1d(i, j)]=1;
	container[_2dto1d(j, i)]=1;
}

void Ensemble_connections::print()
{
	mkdir("params", 0777);

	std::ofstream txtOut;
	txtOut.open("params/"+name+".txt");
	txtOut.precision(8);
	txtOut << N << " ";
	txtOut << k << " ";
	txtOut << proba << " ";
	txtOut << type << " ";
	for (int i = 0; i < container.size(); ++i)
	{
		txtOut << container[i] << " ";
	}
	txtOut.close();
}

void Ensemble_connections::load()
{
	std::ifstream txtIn_check;
	int N_check;
	txtIn_check.open("params.txt");
	txtIn_check >> N_check;
	txtIn_check.close();	
	std::ifstream txtIn;
	txtIn.open("params/"+name+".txt");
	txtIn >> N;
	txtIn >> k;
	txtIn >> proba;
	txtIn >> type;
	container.resize(N*N);
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

std::string Ensemble_connections::get_type()
{
	return type;
}

void Ensemble_connections::print_to_console()
{
	std::cout << name <<": N="<< N << " k="<< k << " proba=" << proba << " " << type<<std::endl;
	
	for (int i = 0; i < N; ++i)
	{
		std::cout << "   ";
		for (int j = 0; j < N; ++j)
		{
			std::cout << container[_2dto1d(i, j)] << " ";
		}

		std::cout << std::endl;
	}
}