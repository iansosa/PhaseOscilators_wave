#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#ifndef ENSEMBLE_CONNECTIONS_H
#define ENSEMBLE_CONNECTIONS_H

class Ensemble_connections
{
	int N;
	std::string name;
	boost::mt19937 &rng;
	std::vector<prec> container;
	int k=2;
	prec proba=1;
	std::string type="global";

	public:

		Ensemble_connections(int in_N, std::string in_name, boost::mt19937 &in_rng, std::string in_type="global") : N(in_N), name(in_name), rng(in_rng), type(in_type), container(N*N) //sets name and all parameters, if all false only first parameter is different from zero
		{
			generate(type); //generates container with first element equal mid
		}

		Ensemble_connections(std::string in_name, boost::mt19937 &in_rng) : name(in_name), rng(in_rng) //sets name and loads parameters from file
		{
			load();
		}

		int size(); //returns size of container

		void resize(int i); //resizes container and generates it again

		void generate(std::string s); //sets elements in container following the instructions in string

		void connect(int i, int j); //connects two elements

		void print(); //prints container in file

		void load(); //loads container from file

		std::string get_type();

		void print_to_console(); //prints container elements in console

		int _2dto1d(int a, int b);

		void _1dto2d(int idx,int &a, int &b);

};

#endif