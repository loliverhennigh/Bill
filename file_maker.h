#include <fstream>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include "Points.h"


using namespace std; 

int seed = 5;

void print_file_rand_low(string file, int particles)
{
	// seed dis shit
	srand(seed);

	// Make da stream
	ofstream writing;
	writing.open(file.c_str());

	// Write it all down
	writing << particles << "\n";

	for (int i = 0; i < 2*particles; i++)
	{
		writing << (float)rand()/((float)RAND_MAX) << " " << (float)rand()/((float)RAND_MAX) << "\n";
	}

	// close dis shit
	writing.close();

}

void print_file(string file, int size, int particles, vector<Point> k)
{
	// Make da stream
	ofstream writing;
	writing.open(file.c_str());

	// Write it all down
	writing << size << "\n" << particles << "\n";

	for (int i = 0; i < particles; i++)
	{
		writing << k[i].x << " " << k[i].y << "\n";
	}

	// close dis shit
	writing.close();
	
}

