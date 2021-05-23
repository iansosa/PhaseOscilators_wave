#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>

#include <boost/numeric/odeint.hpp>
#include <boost/random.hpp>

typedef double prec;

class Dinamic 
{
    prec K;
    int N;

    std::vector<prec> &I;
    std::vector<prec> &F;
    std::vector<std::vector<prec>> &A;
    std::vector<prec> &G;
    std::vector<prec> &W;

	public:
    Dinamic(prec in_K ,int in_N ,std::vector<prec> &in_I ,std::vector<std::vector<prec>> &in_A ,std::vector<prec> &in_F ,std::vector<prec> &in_G ,std::vector<prec> &in_W) : K(in_K) , N(in_N) , I(in_I), A(in_A), F(in_F) , G(in_G), W(in_W)
    { }

    void operator() (const std::vector<prec> &x ,std::vector<prec> &dxdt ,const double t)
    {
    	double sum=0;

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
};

int main()
{
	std::cout << "Hello World" << std::endl;



	return 0;
}