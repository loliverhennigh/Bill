#include <vector>
#include <string>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "bill.h"

using namespace std; 


int main() 
{
	int to = 500;
	int step = 50;
	float dx = .3;
	float pre = .4;

	// generates a random file of balls and prints it to test.txt
	print_file_rand_low("test.txt", 1000000);
	// loads the file just printed
	System a("test.txt");

	// runs the system for time
	for(int i = 0; i < step; i++)
	{
		a.after_time(dx);
	}
	// causes some pertibation
	a.preturb_step(pre);
	// run the system more
	for(int i = 0; i < to-step; i++)
	{
		a.after_time(dx);
	}
	// print all the particle that lie in a small region. This gives us a system with low entropy boundry conditions
	a.System_print_file_low("test1.txt");

	// now run that system and record its entropy for graphing 
	System b("test1.txt");
	for(int i = 0; i < step; i++)
	{
		b.after_time(dx);
	}
	b.preturb_step(pre);
	for(int i = 0; i < to-step; i++)
	{
		b.after_time(dx);
	}
	// print the recorded entropy for graphing
	b.Entropy_print_file("test_pert_50.txt");

	return 0;
}
