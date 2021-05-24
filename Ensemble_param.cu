#include "Ensemble_param.h"


void Ensemble_param::set_rand()
{
	boost::normal_distribution<> gauss(mid,sigma);
    boost::variate_generator< boost::mt19937&, boost::normal_distribution<> > gen(rng,gauss);
    for (int i = 0; i < N; ++i)
    {
    	container[i]=gen();
    }
}

void Ensemble_param::print()
{
	std::ofstream txtOut;
	txtOut.open(name+".txt");
	txtOut.precision(8);
	txtOut << N << " ";
	txtOut << mid << " ";
	txtOut << sigma << " ";
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
	txtIn.open(name+".txt");
	txtIn >> N;
	txtIn >> mid;
	txtIn >> sigma;
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

void Ensemble_param::print_to_console()
{
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