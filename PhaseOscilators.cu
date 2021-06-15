#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>
#include <limits>

#include <sys/stat.h>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/openmp/openmp.hpp>
#include <boost/random.hpp>

typedef double prec;

#include "Ensemble_param.h"
#include "Ensemble_connections.h"
#include "Ensemble_Dinamic.h"
#include "Ensemble_evolve.h"


void calc_U(Evolve &e, Dinamic &P, prec dx, std::vector<prec> &U)
{
	for (int i = 1; i < U.size(); ++i)
	{
		prec aux=e.Amp(i,P);
		if(aux>0)
		{
			U[i]=log(aux);
		}
	}

	for (int i = 1; i < U.size()-1; ++i)
	{
		U[i]=(U[i+1]-U[i])/dx;
	}
}

void printU(Evolve &e, Dinamic &P, prec dx)
{
	std::vector<prec> U(e.size(),0);

	calc_U(e,P,dx,U);

	std::string name="U";
	mkdir("out", 0777);
	std::ofstream txtOut;
	txtOut.open("out/"+name+".txt");
	txtOut.precision(std::numeric_limits< prec >::max_digits10);

	for (int i = 1; i < U.size()-1; ++i)
	{
			txtOut << i << " " << U[i] << std::endl;
	}
	txtOut.close();	
}

void print_U_stats(std::vector<std::vector<prec>> &U)
{
	std::vector<prec> U_med(U[0].size(),0);
	std::vector<prec> U_sigm(U[0].size(),0);
	for (int i = 0; i < U[0].size(); ++i) //N
	{
		int norma=0;
		for (int j = 0; j < U.size(); ++j) //promedio
		{
			if(U[j][i]>-1000 && U[j][i]<500)
			{
				norma++;
				U_med[i]=U_med[i]+U[j][i];
			}	
		}
		if(norma>0)
		{
			U_med[i]=U_med[i]/norma;
		}
	}

	for (int i = 0; i < U[0].size(); ++i) //N
	{
		int norma=0;
		for (int j = 0; j < U.size(); ++j) //promedio
		{
			if(U[j][i]>-1000 && U[j][i]<500)
			{
				norma++;
				U_sigm[i]=pow(U[j][i]-U_med[i],2);
			}	
		}
		if(norma>0)
		{
			U_sigm[i]=U_sigm[i]/norma;
		}
		U_sigm[i]=sqrt(U_sigm[i]);
	}
	std::string name="U_stats";
	mkdir("out", 0777);
	std::ofstream txtOut;
	txtOut.open("out/"+name+".txt");
	txtOut.precision(std::numeric_limits< prec >::max_digits10);

	for (int i = 0; i < U_med.size(); ++i)
	{
			txtOut << i << " " << U_med[i] << " " << U_sigm[i] << std::endl;
	}
	txtOut.close();	
}

void promedio_U()
{
 	boost::mt19937 rng(static_cast<unsigned int>(std::time(0)));

 	int N;
 	std::cout << "N: ";
 	std::cin >> N;

  	int prom;
 	std::cout << "Prom: ";
 	std::cin >> prom;

 	std::vector<std::vector<prec> > U(prom);

    int chunk_size = N/omp_get_max_threads();
    omp_set_schedule( omp_sched_static , chunk_size );

    prec conv_factor=5;
    prec dx=0.01;
    dx=dx/conv_factor;

    prec K=conv_factor*conv_factor*N;
    prec I=0.1;
    prec G=5;
	Dinamic P(N ,K ,rng ,I ,I*0.01 ,1 ,0 ,G ,G*0.01 ,1 ,0,"chain"); //N, K, rng, I, sigm_I, F, sigm_F, G, sigm_G, W, sigm_W
	P.init_I_type_rand(true,false);

	Evolve e(N, rng);
	for (int i = 0; i < prom; ++i)
	{
		std::cout << "Iter: (" << i+1 << "/" << prom << ")" << std::endl;
		U[i].resize(N);
		e.run(P,0,50000,100);
		e.calc_convergence();
		calc_U(e,P,dx,U[i]);
		e.clean();
		P.print_params();

		P.generate();
	}
	e.print();

	print_U_stats(U);
}

int main()
{
 	boost::mt19937 rng(static_cast<unsigned int>(std::time(0)));

 	int N;
 	std::cout << "N: ";
 	std::cin >> N;
 	int chunk_size = N/omp_get_max_threads();
 	omp_set_schedule( omp_sched_static , chunk_size );

    prec I=0.1;
    prec G=0.25;
    prec K=N;
	Dinamic P(N ,K ,rng ,I ,I*0 ,1000 ,0 ,G ,G*0 ,1 ,0,"chain"); //N, K, rng, I, sigm_I, F, sigm_F, G, sigm_G, W, sigm_W
	Evolve e(N, rng);

	e.run(P,0,200000,1000);
	e.print();
	e.calc_convergence();
	printU(e,P,0.01);

	return 0;
}