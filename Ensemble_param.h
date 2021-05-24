#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#ifndef ENSEMBLE_PARAM_H
#define ENSEMBLE_PARAM_H

#include <boost/numeric/odeint.hpp>
#include <boost/random.hpp>

typedef double prec;

class Ensemble_param
{
	int N;
	std::string name;
	boost::mt19937 &rng;

	public:
		std::vector<prec> container;
		prec mid;
		prec sigma;
		Ensemble_param(int in_N, std::string in_name, prec in_mid, prec in_sigma, boost::mt19937 &in_rng) : N(in_N), name(in_name), rng(in_rng), mid(in_mid), sigma(in_sigma), container(N)
		{
			set_rand();
		}

		Ensemble_param(std::string in_name, boost::mt19937 &in_rng) : name(in_name), rng(in_rng)
		{
			load();
		}

		prec& operator[](int i);

		void set_rand();

		void print();

		void load();

		void print_to_console();
};

#endif