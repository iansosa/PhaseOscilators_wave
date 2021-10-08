#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#ifndef ENSEMBLE_CONNECTIONS_H
#define ENSEMBLE_CONNECTIONS_H

class Ensemble_connections
{
	int N; //numbrt of connections
	std::string name; //internal name of variable
	boost::mt19937 &rng; //random variable generator
	std::vector<prec> container; //conection container
	int k=2; //connection grade of ensemble, not used
	prec proba=1; //probability of activation of conection, not used
	std::string type="global"; //type of ensemble

	void set_type_id(); //set the type_id according to type name

	public:

		int type_id=-1; //numerical value representing of type

		Ensemble_connections(int in_N, std::string in_name, boost::mt19937 &in_rng, std::string in_type="global") : N(in_N), name(in_name), rng(in_rng), type(in_type), container(N*N) //sets name and all parameters, if all false only first parameter is different from zero
		{
			generate(type); //generates container with first element equal mid
			set_type_id();
		}

		Ensemble_connections(std::string in_name, boost::mt19937 &in_rng) : name(in_name), rng(in_rng) //sets name and loads parameters from file
		{
			load();
			set_type_id();
		}

		prec operator() (int a, int b)
		{
			return container[_2dto1d(a,b)];
		}

		int size(); //returns size of container
		void resize(int i); //resizes container and generates it again
		void generate(std::string s); //sets elements in container following the instructions in string
		void connect(int i, int j); //connects two elements
		void print(); //prints container in file
		void load(); //loads container from file
		std::string get_type(); //returns type
		void print_to_console(); //prints container elements in console
		bool TypeIs(std::string s); //true if type contains the word s
		int _2dto1d(int a, int b); //transform a and b to a 1D coordinate
		void _1dto2d(int idx,int &a, int &b); //using idx sets a and b to the corresponding 2D coordinates
};

#endif