#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#ifndef ENSEMBLE_PARAM_H
#define ENSEMBLE_PARAM_H

class Ensemble_param
{
	int N;
	std::string name;
	boost::mt19937 &rng;
	std::vector<prec> container;
	prec mid;
	prec sigma;
	std::string type="gauss";

	public:

		Ensemble_param(int in_N, std::string in_name, prec in_mid, boost::mt19937 &in_rng, bool all=true) : N(in_N), name(in_name), rng(in_rng), mid(in_mid), container(N) //sets name and all parameters, nothing random, if all false only first parameter is different from zero
		{
			sigma=0;
			type_rand(all);
			generate(false); //generates container with first element equal mid
		}

		Ensemble_param(int in_N, std::string in_name, prec in_mid, prec in_sigma, boost::mt19937 &in_rng, bool all=true) : N(in_N), name(in_name), rng(in_rng), mid(in_mid), sigma(in_sigma), container(N) //sets name and all parameters, if all false only first parameter is different from zero
		{
			type_rand(all);
			generate(false); //generates container with first element equal mid
		}

		Ensemble_param(std::string in_name, boost::mt19937 &in_rng) : name(in_name), rng(in_rng) //sets name and loads parameters from file
		{
			load();
		}

		void dist(prec in_mid, prec in_sigma=0); //sets medium value and sigma

		int size(); //returns size of container

		void resize(int i); //resizes container

		prec& operator[](int i); //returns element i from container

		void generate(bool first=false); //sets elements in container with a gaussian distribution(mid,sigma). If argument true then first element is also random

		void print(); //prints container in file

		void load(); //loads container from file

		void type_rand(bool=true); //if true set type to gauss, if false type is set to ofirst

		void print_to_console(); //prints container elements in console


};

#endif