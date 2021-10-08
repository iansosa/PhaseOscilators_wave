#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>
#include <limits>
#include <cmath>

#include <sys/stat.h>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/openmp/openmp.hpp>
#include <boost/random.hpp>

typedef long double prec;

#include "../../Ensemble_param.h"
#include "../../Ensemble_connections.h"
#include "../../Ensemble_Dinamic.h"
#include "../../Ensemble_evolve.h"

int period(std::vector<int> &v)
{
	bool equal=true;
	for (int i = 1; i < v.size()/2; ++i)
	{
		equal=true;
		for (int j = 1; j < v.size()-i; ++j)
		{
			if(v[j]!=v[j+i])
			{
				equal=false;
				break;
			}
		}
		if(equal==true)
		{
			return i;
		}
	}
	return -1;
} 

void reevaluate_F(int N,Dinamic &P,prec F_start, prec dF)
{
	for (int i = 0; i < N; ++i)
	{
		P.new_va_F(i,F_start+i*dF);			
	}
}

void reevaluate_F_lyapunov(int N,Dinamic &P,prec F_start, prec dF)
{
	for (int i = 0; i < N; ++i)
	{
		if(i%2==0)
		{
			P.new_va_F(i,F_start+i*dF);	
		}
		else
		{
			P.new_va_F(i,F_start+(i-1)*dF);	
		}
	}
}

void calc_save_n(prec G,std::vector <int> &n, Evolve &e, prec F_start, prec dF, prec t_save, prec dt, bool secondary=false)
{
	for (int i = 0; i < n.size(); ++i)
	{
		n[i]=(int)((e.x(i,t_save/dt)-e.x(i,0))/(2.0*M_PI));
		std::cout << i << " F=" << F_start+i*dF << " n=" << n[i] << " x=" << e.x(i,t_save/dt)-e.x(i,0) <<std::endl;
	}

	std::string G_s=std::to_string(G);
	G_s = G_s.substr(0, 6); 
	std::string filename="out/n_"+G_s+".txt";
	if(secondary==true)
	{
		filename="out/n_"+G_s+"_b.txt";
	}
	std::ofstream txtOut;
	txtOut.precision(std::numeric_limits< prec >::max_digits10);
	txtOut.open(filename);
	for (int j = 0; j < n.size(); ++j)
	{
		txtOut << j << " " << F_start+j*dF << " " << n[j] << std::endl;
	}
	txtOut.close();
}


void calc_save_b(prec G,std::vector<std::vector<int>>&b, Evolve &e, prec F_start, prec dF, prec t_save, prec dt, bool secondary=false)
{
	for (int i = 0; i < b.size(); ++i)
	{
		int prev=-1;
		for (int j = 0; j < (int) (t_save/dt); ++j)
		{
			int curr=j*dt/(2.0*M_PI);
			if(curr!=prev)
			{
				prev=curr;
				int bites=(fabs(e.x(i,j)-M_PI)/(2.0*M_PI));
				b[i].push_back(bites);
			}
		}
	}

	for (int i = 0; i < b.size(); ++i)
	{
		for (int j = b[i].size()-1; j >= 1; --j)
		{
			b[i][j]=b[i][j]-b[i][j-1];
		}
	}
	std::string G_s=std::to_string(G);
	G_s = G_s.substr(0, 6); 
	std::string filename="out/b_"+G_s+".txt";
	if(secondary==true)
	{
		filename="out/b_"+G_s+"_b.txt";
	}
	std::ofstream txtOut;
	txtOut.precision(std::numeric_limits< prec >::max_digits10);
	txtOut.open(filename);
	for (int i = 0; i < b.size(); ++i)
	{
		txtOut << F_start+i*dF << " ";
		for (int j = 1; j < b[i].size(); ++j)
		{
			txtOut << b[i][j]<< " ";
		}
		txtOut << std::endl;
	}
	txtOut.close();
}

void calc_save_map(prec G,int N, Evolve &e, prec F_start, prec dF, prec t_save, prec dt, bool secondary=false)
{
	std::vector<std::vector<prec>> map_x(N);
	std::vector<std::vector<prec>> map_v(N);
	for (int i = 0; i < N; ++i)
	{
		int prev=-1;
		for (int j = 0; j < (int) (t_save/dt); ++j)
		{
			int curr=(j*dt+1.0)/(2.0*M_PI);
			if(curr!=prev)
			{
				prev=curr;
				int n= e.x(i,j)/(2.0*M_PI);
				map_x[i].push_back(e.x(i,j)-n*2.0*M_PI);
				map_v[i].push_back(e.v(i,j));
			}
		}
	}

	std::string G_s=std::to_string(G);
	G_s = G_s.substr(0, 6); 
	std::string filename="out/map_"+G_s+".txt";
	if(secondary==true)
	{
		filename="out/map_"+G_s+"_b.txt";
	}
	std::ofstream txtOut;
	txtOut.precision(std::numeric_limits< prec >::max_digits10);
	txtOut.open(filename);
	for (int i = 0; i < map_x.size(); ++i)
	{
		txtOut << F_start+i*dF << " ";
		for (int j = 1; j < map_x[i].size(); ++j)
		{
			txtOut << map_x[i][j]<< " ";
		}
		txtOut << std::endl;

		txtOut << F_start+i*dF << " ";
		for (int j = 1; j < map_v[i].size(); ++j)
		{
			txtOut << map_v[i][j]<< " ";
		}
		txtOut << std::endl;
	}
	txtOut.close();

	filename="out/per_"+G_s+".txt";
	txtOut.open(filename);
	for (int i = 0; i < map_v.size(); ++i)
	{
		for (int j = 1; j < map_v[i].size(); ++j)
		{
			txtOut << F_start+i*dF << " " << map_v[i][j]<< std::endl;
		}
	}
	txtOut.close();
}


void calc_save_p(prec G,std::vector<std::vector<int>>&b, Evolve &e, prec F_start, prec dF, prec t_save, prec dt, bool secondary=false)
{
	std::string G_s=std::to_string(G);
	G_s = G_s.substr(0, 6); 
	std::string filename="out/p_"+G_s+".txt";
	if(secondary==true)
	{
		filename="out/p_"+G_s+"_b.txt";
	}
	std::ofstream txtOut;
	txtOut.precision(std::numeric_limits< prec >::max_digits10);
	txtOut.open(filename);
	for (int i = 0; i < b.size(); ++i)
	{
		txtOut << F_start+i*dF << " " << period(b[i]) << std::endl;
	}
	txtOut.close();
}

void set_init_lyapunov(boost::mt19937 &rng,Evolve &e,prec t_save, prec dt)
{
    boost::uniform_real<> unif( 0, 2*M_PI );
    boost::variate_generator< boost::mt19937&, boost::uniform_real<> > gen( rng , unif );

	std::vector<prec> x_in(e.size());
	std::vector<prec> v_in(e.size());

	for (int i = 0; i < e.size(); ++i)
	{
		prec phase=gen();
		if(i%2==0)
		{
			x_in[i]=e.x(i,t_save/dt);
			v_in[i]=e.v(i,t_save/dt);
		}
		else
		{
			x_in[i]=e.x(i-1,t_save/dt)+0.0000001*cos(phase);
			v_in[i]=e.v(i-1,t_save/dt)+0.0000001*sin(phase);
		}
	}
    e.set_init(x_in,v_in);
}

void calc_lyapunov(std::vector<prec> &l,std::vector<int> &l_prom,Evolve &e,prec t_save, prec dt,prec T)
{
	for (int i = 0; i < l.size(); ++i)
	{
		prec x1=e.x(2*i,T/dt);
		prec x2=e.x(2*i+1,T/dt);
		prec v1=e.v(2*i,T/dt);
		prec v2=e.v(2*i+1,T/dt);
		prec dx=fabs(x1-x2);
		int n=(int)(dx/(2.0*M_PI));
		dx=dx-n*2*M_PI;
		if(dx>=M_PI)
		{
			dx=2*M_PI-dx;
		}
		prec dv=v1-v2;
		prec lya=(log(sqrt(pow(dx*10000000,2)+pow(dv*10000000,2)))-log(0.0000001*10000000))/(T);
		if(std::isinf(lya)==false)
		{
			l[i]=l[i]+lya;
			l_prom[i]=l_prom[i]+1;
		}

	}

}

void measure_lyapunov(int N,prec F_start, prec F_end, prec G, prec K, prec evolution_time, int traject_num)
{
 	boost::mt19937 rng(static_cast<unsigned int>(std::time(0)));

    boost::uniform_real<> unif( 0, 2*M_PI );
    boost::variate_generator< boost::mt19937&, boost::uniform_real<> > gen( rng , unif );

 	N=2*N;

 	prec dt=0.01;
 	prec t_start=0;
 	prec t_end=10000;
 	prec t_save=evolution_time;

    prec I=1;

	Dinamic P(N ,K ,rng ,I ,I*0 , 1.3794,0 ,G ,G*0 ,1,0,"solid_chain"); //N, K, rng, I, sigm_I, F, sigm_F, G, sigm_G, W, sigm_W
	Evolve e(N, rng);


	prec dF=(F_end-F_start)/N;
	reevaluate_F_lyapunov(N,P,F_start,dF);
	P.print_params();
	e.run(P,t_start,t_end,t_save,dt);
	
	t_start=t_end;
	t_end=t_end+t_save;
	set_init_lyapunov(rng,e,t_save,dt);
	e.clean();

	e.print_init();
	std::vector<long double> l(e.size()/2,0);
	std::vector<int> l_prom(e.size()/2,0);

	for (int i = 0; i < traject_num; ++i)
	{
		e.run(P,t_start,t_end,t_save,dt);
		calc_lyapunov(l,l_prom,e,t_save,dt,evolution_time);

		t_start=t_end;
		t_end=t_end+t_save+gen();
		set_init_lyapunov(rng,e,t_save,dt);
		if(i!= traject_num-1)
		{
			e.clean();
		}
	}

	std::string G_s=std::to_string(G);
	G_s = G_s.substr(0, 6); 
	std::string filename="out/lyapunov_"+G_s+".txt";

	std::ofstream txtOut;
	txtOut.precision(std::numeric_limits< prec >::max_digits10);
	txtOut.open(filename);
	for (int i = 0; i < l.size(); ++i)
	{
		txtOut << F_start+2*i*dF << " " <<l[i]/l_prom[i]<< std::endl;
	}
	txtOut.close();

	filename="out/lyapunov_d_"+G_s+".txt";
	txtOut.open(filename);
	prec max=0;
	for (int i = 0; i < l.size(); ++i)
	{
		if(l[i]>=max)
		{
			max=l[i];
		}
	}
	for (int i = 0; i < l.size(); ++i)
	{
		if(l[i]>max/5)
		{
			txtOut << F_start+2*i*dF << " " << 1<< std::endl;
		}
		else
		{
			txtOut << F_start+2*i*dF << " " << 0<< std::endl;
		}
	}
	txtOut.close();

	//calc_save_n(G,n,e,F_start,dF,t_save,dt,false);
	//calc_save_b(G,b,e,F_start,dF,t_save,dt,false);
	//calc_save_p(G,b,e,F_start,dF,t_save,dt,false);
	//calc_save_map(G,N,e,F_start,dF,t_save,dt,false);

	e.print(false);
	e.clean();
}

void measure_chaos(int N,prec F_start, prec F_end, prec G, prec K)
{
 	boost::mt19937 rng(static_cast<unsigned int>(std::time(0)));

 	prec dt=0.01;
 	prec t_start=0;
 	prec t_end=20000;
 	prec t_save=10000;

    prec I=1;

	Dinamic P(N ,K ,rng ,I ,I*0 , 1.3794,0 ,G ,G*0 ,1,0,"solid_chain"); //N, K, rng, I, sigm_I, F, sigm_F, G, sigm_G, W, sigm_W
	Evolve e(N, rng);

	std::vector <int> n(N,0);
	std::vector <std::vector <int> > b(N);
	prec dF=(F_end-F_start)/N;
	reevaluate_F(N,P,F_start,dF);
	P.print_params();
	e.run(P,t_start,t_end,t_save,dt);

	calc_save_n(G,n,e,F_start,dF,t_save,dt,false);
	calc_save_b(G,b,e,F_start,dF,t_save,dt,false);
	calc_save_p(G,b,e,F_start,dF,t_save,dt,false);
	calc_save_map(G,N,e,F_start,dF,t_save,dt,false);

	//natural_scale(10,0.01);
	e.print(false);
	e.clean();
}

prec FofP(prec G,prec A0,prec A1,prec A2,prec A3,prec A4)
{
	return A0+A1*G+A2*G*G+A3*G*G*G+A4*G*G*G*G;
}

void param_iter(int N,int N_P,prec P_start, prec P_end)
{
	prec F_start;
	prec F_end;
 	prec dP=(P_end-P_start)/N_P;

 	std::ofstream txtOut;
 	for (int i = 117; i <= N_P; ++i)
 	{
 		F_start=FofP(P_start+dP*i,0,0.4,0,0.0,0.0);
 		F_end=FofP(P_start+dP*i,2.5,0.7,0.0,0.0,0.0);
 		std::cout << "iter: " << i << "/" << N_P << std::endl;
 		std::cout << "G=" << P_start+dP*i << std::endl;
 		std::cout << "F_start=" << F_start << std::endl;
 		std::cout << "F_end=" << F_end << std::endl;
 		//measure_chaos(N,F_start,F_end,P_start+dP*i,1);
 		measure_lyapunov(N,F_start,F_end,P_start+dP*i,1,50.0,100);
		
		txtOut.open("checkpoint.txt");
		txtOut << i  << std::endl;
		txtOut.close();
 	}
}


int main()
{
 	int N;
 	std::cout << "N: ";
 	std::cin >> N;

 	int nProcessors=omp_get_max_threads();
  	nProcessors=nProcessors-2;
  	//nProcessors=4;
    int chunk_size = N/nProcessors;
    omp_set_num_threads(nProcessors);
    omp_set_schedule( omp_sched_static , chunk_size );

 	prec F_start=1.7;
 	prec F_end=1.9;
 	prec G=1.5;
 	prec K=1;

 	//measure_chaos(N,F_start,F_end,G,K);
 	//param_iter(N,500,0.1,5);
 	measure_lyapunov(N,F_start,F_end,G,K,251.32,500);

	return 0;
}