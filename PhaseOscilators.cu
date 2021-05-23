#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>

#include <boost/numeric/odeint.hpp>
#include <boost/random.hpp>

typedef double prec;

class Ensemble_param
{
	int N;
	std::string name;
	boost::mt19937 *rng;

	public:
		std::vector<prec> container;
		prec mid;
		prec sigma;
		Ensemble_param(int in_N, std::string in_name, prec in_mid, prec in_sigma, boost::mt19937 &in_rng) : N(in_N), name(in_name), rng(&in_rng), mid(in_mid), sigma(in_sigma), container(N)
		{
			set();
		}

		Ensemble_param(int in_N, std::string in_name) : N(in_N), name(in_name), container(N), rng(NULL)
		{
			load();
		}

		prec& operator[](int i)
		{
			return container[i];
		}

		void set()
		{
			boost::normal_distribution<> gauss(mid,sigma);
    		boost::variate_generator< boost::mt19937&, boost::normal_distribution<> > gen(*rng,gauss);

    		for (int i = 0; i < N; ++i)
    		{
    			container[i]=gen();
    		}
		}

		void print()
		{
			std::ofstream txtOut;
			txtOut.open(name+".txt");
			txtOut.precision(8);
			for (int i = 0; i < N; ++i)
			{
				txtOut << container[i] << " ";
			}
			txtOut.close();
		}

		void load()
		{
			std::ifstream txtIn;
			txtIn.open(name+".txt");
			for (int i = 0; i < N; ++i)
			{
				txtIn >> container[i];
			}
			txtIn.close();
		}

		void print_to_console()
		{
			for (int i = 0; i < N; ++i)
			{
				std::cout << container[i] << " ";
			}
			std::cout << std::endl;
		}

};

class Dinamic 
{
    prec K;

    Ensemble_param I;
    Ensemble_param F;
    std::vector<std::vector<prec>> A;
    Ensemble_param G;
    Ensemble_param W;

	public:
		int N;

    	Dinamic(int in_N , prec in_K) : N(in_N), K(in_K), I(in_N,"I"), F(in_N,"F"), A(in_N), G(in_N,"G"), W(in_N,"W")
    	{ }

    	Dinamic(int in_N , prec in_K, boost::mt19937 &in_rng, prec mid_I, prec sigma_I, prec mid_F, prec sigma_F, prec mid_G, prec sigma_G, prec mid_W, prec sigma_W) : N(in_N), K(in_K), I(in_N,"I",mid_I,sigma_I,in_rng), F(in_N,"F",mid_F,sigma_F,in_rng), A(in_N), G(in_N,"G",mid_G,sigma_G,in_rng), W(in_N,"W",mid_W,sigma_W,in_rng)
    	{ }

    	void operator() (const std::vector<prec> &x ,std::vector<prec> &dxdt ,const double t)
    	{
    	    for (int i = 0; i < N; ++i)
    	    {
    	    	prec sum=interaction_sum_chain(i,x);

        		dxdt[2*i]=x[2*i+1];
				dxdt[2*i+1]=sum/I[i]+F[i]*sin(W[i]*t-x[2*i])/I[i]-(G[i]/I[i])*x[2*i+1];      	
        	}
    	}

    	prec interaction_sum_all(int id, const std::vector<prec> &x)
    	{
    		prec sum=0;
    		for (int i = 0; i < N; ++i)
    		{
    			sum=sum+A[id][i]*sin(x[2*i]-x[2*id])/N;
    		}
    		return sum*K;
    	}

    	prec interaction_sum_chain(int id,const std::vector<prec> &x)
    	{
    		prec sum=0;
    	    if(id==0)
    	    {
    	    	sum=sum+A[id][id+1]*sin(x[2*(id+1)]-x[2*id])/N;
    	    }
    	    if(id==N-1)
    	    {
    	    	sum=sum+A[id][id-1]*sin(x[2*(id-1)]-x[2*id])/N;
    	    }
    	    if(id>0 && id<N-1)
    	    {
    	    	sum=sum+A[id][id+1]*sin(x[2*(id+1)]-x[2*id])/N;
    	    	sum=sum+A[id][id-1]*sin(x[2*(id-1)]-x[2*id])/N;
    	    }
    	   
    		return sum*K;
    	}

    	void print_params()
    	{
    		I.print();
    		F.print();
    		G.print();
    		W.print();
    	}

    	void print_params_to_console()
    	{
    		I.print_to_console();
    		F.print_to_console();
    		G.print_to_console();
    		W.print_to_console();
    	}
};

int main()
{
 	boost::mt19937 rng(static_cast<unsigned int>(std::time(0)));

	std::cout << "Hello World" << std::endl;


	Dinamic P(10, 1);
	P.print_params_to_console();

	return 0;
}