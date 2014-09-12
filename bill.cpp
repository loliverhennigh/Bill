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

	print_file_rand_low("test.txt", 1000000);
	System a("test.txt");
//	printf("%g \n", a.entropy);

	for(int i = 0; i < step; i++)
	{
		a.after_time(dx);
	}
	a.preturb_step(pre);
	for(int i = 0; i < to-step; i++)
	{
		a.after_time(dx);
	}
	a.System_print_file_low("test1.txt");


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

	b.Entropy_print_file("test_pert_50.txt");

	return 0;
}
