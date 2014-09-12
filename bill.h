
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "Points.h"

using namespace std; 

// seed for random stuff
int seed = 31;




// print functions
void print_matrix(gsl_matrix * k, int size);
void print_vector(gsl_vector * v, int size);
void print_file_rand_low(string file, int size, int particles);
void print_file(string file, int size, vector<Point> k);
void print_file_vector(string file, vector<float> p);
void print_point(Point p);
void print_vector_vector(vector<vector<int> > p);
void print_file_vector(string file, vector<float> k);

// the dusy
Point after_time_single( Point p, float t);


class System {
	public:
	// number of particles
	int particles;
	// entropy
	float entropy;
	// vector of points
	vector<Point> points;
	// vector of starting points
	vector<Point> start_points;
	// stores the entorpy
	vector<float> entropys;
	// middle bit
	vector<vector<int> > middle_box;
	vector<vector<int> > left_circle;
	vector<vector<int> > right_circle;



	// constroctors read in the matrix	(uses cin haha sorry) 
	System(string file);
	// print it for testing
	void System_print();
	// print points for testing
	void Points_print();
	// print ensemble for testing
	void Ensemble_print();
	// calc new ensemble
	void renew_ensemble();
	// find distibution at time t
	void after_time(float t);
	// preturb it!
	void preturb_step(float dp);
	// calc entropy
	void calc_entropy();
	// print to file
	void System_print_file(string file);
	// print low to file
	void System_print_file_low(string file);
	// print entropy vector
	void Entropy_print_file(string file);
	// one time step
	void one_step_coupling();

};


// constructor guy
System::System(string file)
{
	// helpful later
	float store_x;
	float store_y;
	float store_x_v;
	float store_y_v;

	// init the box counters
	for (int i = 0; i < 10; i++)
	{
		vector<int> temp;
		for (int j = 0; j < 20; j++)
		{
			temp.push_back(0);
		}

		middle_box.push_back(temp);
	}
	
	for (int i = 0; i < 5; i++)
	{
		vector<int> a;
		for (int k = 0; k < 28; k++)
		{
			a.push_back(0);
		}
		left_circle.push_back(a);
		right_circle.push_back(a);
	}
	
	// creat file stream
	ifstream reading;
	reading.open(file.c_str());

	// just in case allocations
	entropy = 0;
	particles = 0;

	// start the porly coded madness!
	reading >> particles;
	
	for (int i = 0; i < particles; i++)
	{
		reading >> store_x >> store_y >> store_x_v >> store_y_v;
		Point p = {store_x, store_y, store_x_v, store_y_v};
		points.push_back(p);
		start_points.push_back(p);
	}
		
	// remember to close!
	reading.close();

	
	renew_ensemble();
	//print_vector_vector(middle_box);
	//print_vector_vector(left_circle);
	//print_vector_vector(right_circle);
	calc_entropy();
	entropys.push_back(entropy);

}
// print guy
void System::System_print()
{
	printf ("particles \n");
	Points_print();
	printf ("entropy %g \n", entropy);
}
// print paritcles
void System::Points_print()
{
	for (int i = 0; i < particles; i++)
	{
		printf("(%g, %g) velocity (%g, %g) \n", points[i].x, points[i].y, points[i].v_x, points[i].v_y);
	}
}
// print maticx
void System::Ensemble_print()
{
}

// new ensemble guy
void System::renew_ensemble()
{
	// make zero, bad
	for(int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 20; j++)
		{
			middle_box[i][j] = 0;
		}
	}
	
	for (int i = 0; i < 5; i++)
	{
		for (int k = 0; k < 3; k++)
		{
			left_circle[i][k] = 0;
			right_circle[i][k] = 0;
		}
		for (int k = 0; k < 9; k++)
		{
			left_circle[i][k] = 0;
			right_circle[i][k] = 0;
		}
		for (int k = 0; k < 16; k++)
		{
			left_circle[i][k] = 0;
			right_circle[i][k] = 0;
		}
		for (int k = 0; k < 22; k++)
		{
			left_circle[i][k] = 0;
			right_circle[i][k] = 0;
		}
		for (int k = 0; k < 28; k++)
		{
			left_circle[i][k] = 0;
			right_circle[i][k] = 0;
		}

	}

	for (int i = 0; i < particles; i++)
	{
		if (points[i].x >= -1 && points[i].x <= 1 && points[i].y >= -1 && points[i].y <= 1)
		{
			int spot_x_mid = int((points[i].x+1)/.2);
			int spot_y_mid = int((points[i].y+1)/.1);
			if(spot_x_mid > 9) {spot_x_mid = 9;}
			if(spot_x_mid < 0) {spot_x_mid = 0;}
			if(spot_y_mid > 19) {spot_y_mid = 19;}
			if(spot_y_mid < 0) {spot_y_mid = 0;}
			middle_box[spot_x_mid][spot_y_mid]++;
		}
		if (points[i].x < -1)
		{
			// 28 layer 1 radious > .8
			// 22 layer 2 radious > .6 
			// 16 layer 3 radious > .4
			// 9 layer 4 radious > .2
			// 3 layer 5 radious > 0
			float r = sqrt((points[i].x+1)*(points[i].x+1)+points[i].y*points[i].y);
			float theta = atan(points[i].y/(points[i].x+1)) + 3.14159/2;
			int spot_r_left = int(r/.2);
			if(spot_r_left > 4) {spot_r_left = 4;}
			int spot_theta_left = 0;	
			if (spot_r_left == 0)
			{	
				spot_theta_left = int(theta/1.0471);
				if(spot_theta_left > 2){ spot_theta_left = 2;}
			}
			if (spot_r_left == 1)
			{	
				spot_theta_left = int(theta/.3490655);
				if(spot_theta_left > 8){ spot_theta_left = 8;}
			}
			if (spot_r_left == 2)
			{	
				spot_theta_left = int(theta/.19634);
				if(spot_theta_left > 14){ spot_theta_left = 15;}
			}
			if (spot_r_left == 3)
			{	
				spot_theta_left = int(theta/.142799);
				if(spot_theta_left > 21){ spot_theta_left = 21;}
			}
			if (spot_r_left == 4)
			{	
				spot_theta_left = int(theta/.112199);
				if(spot_theta_left == 28){ spot_theta_left = 27;}
			}
			if(spot_theta_left < 0){ spot_theta_left = 0;}

			left_circle[spot_r_left][spot_theta_left]++;
		}
		if (points[i].x > 1)
		{
			// 28 layer 1 radious > .8
			// 22 layer 2 radious > .6 
			// 16 layer 3 radious > .4
			// 9 layer 4 radious > .skf2
			// 3 layer 5 radious > 0
			float r = sqrt((points[i].x-1)*(points[i].x-1)+points[i].y*points[i].y);
			float theta = atan(points[i].y/(points[i].x-1)) + 3.14159/2;
			int spot_r_right = int(r/.2);
			if(spot_r_right > 4) {spot_r_right = 4;}
			int spot_theta_right = 0;
			if (spot_r_right == 0)
			{	
				spot_theta_right = int(theta/1.0471);
				if(spot_theta_right == 3){ spot_theta_right = 2;}
			}
			if (spot_r_right == 1)
			{	
				spot_theta_right = int(theta/.3490655);
				if(spot_theta_right == 9){ spot_theta_right = 8;}
			}
			if (spot_r_right == 2)
			{	
				spot_theta_right = int(theta/.19634);
				if(spot_theta_right == 16){ spot_theta_right = 15;}
			}
			if (spot_r_right == 3)
			{	
				spot_theta_right = int(theta/.142799);
				if(spot_theta_right == 22){ spot_theta_right = 21;}
			}
			if (spot_r_right == 4)
			{	
				spot_theta_right = int(theta/.112199);
				if(spot_theta_right == 28){ spot_theta_right = 27;}
			}
			if(spot_theta_right < 0){ spot_theta_right = 0;}
			right_circle[spot_r_right][spot_theta_right]++;
		}
	}	

}

// take a time step
void System::after_time(float t)
{
	for (int i = 0; i < particles; i++)
	{
		points[i] = after_time_single(points[i], t);
	}

	renew_ensemble();
	calc_entropy();
	entropys.push_back(entropy);
	//print_vector_vector(middle_box);
	//print_vector_vector(left_circle);
	//print_vector_vector(right_circle);
}

// the preturb step!
void System::preturb_step(float dp)
{
	for (int i = 0; i < particles; i++)
	{
		points[i].v_x = points[i].v_x + dp;

		//float n_x = 4181*(points[i].x+1)+6765*(points[i].y+1);
		//float n_y = 10946*(points[i].y+1)+6765*(points[i].x+1);
		//n_x = n_x - 2*int(n_x/2.0);
		//n_y = n_y - 2*int(n_y/2.0);
		//if(n_x > 0 && n_x < 2 && n_y > 0 && n_y < 2)
		//{
		//	points[i].x = n_x-1;
		//	points[i].y = n_y-1;
		//}

	}

}

void System::calc_entropy()
{
	entropy = 0.0;
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 10; j++)
		{
			float store_p = float(middle_box[i][j])/float(particles);
			if(store_p > 0)
			{
				entropy = entropy - store_p * log(store_p);
			}
		}
	}
	
	for (int i = 0; i < 5; i++)
	{
		for (int k = 0; k < 3; k++)
		{
			float store_l = float(left_circle[i][k])/float(particles);
			float store_r = float(right_circle[i][k])/float(particles);
			if(store_l > 0)
			{
				entropy = entropy - store_l * log(store_l);
			}
			if(store_r > 0)
			{
				entropy = entropy - store_r * log(store_r);
			}
		}
		for (int k = 0; k < 9; k++)
		{
			float store_l = float(left_circle[i][k])/float(particles);
			float store_r = float(right_circle[i][k])/float(particles);
			if(store_l > 0)
			{
				entropy = entropy - store_l * log(store_l);
			}
			if(store_r > 0)
			{
				entropy = entropy - store_r * log(store_r);
			}
		}
		for (int k = 0; k < 16; k++)
		{
			float store_l = float(left_circle[i][k])/float(particles);
			float store_r = float(right_circle[i][k])/float(particles);
			if(store_l > 0)
			{
				entropy = entropy - store_l * log(store_l);
			}
			if(store_r > 0)
			{
				entropy = entropy - store_r * log(store_r);
			}
		}
		for (int k = 0; k < 22; k++)
		{
			float store_l = float(left_circle[i][k])/float(particles);
			float store_r = float(right_circle[i][k])/float(particles);
			if(store_l > 0)
			{
				entropy = entropy - store_l * log(store_l);
			}
			if(store_r > 0)
			{
				entropy = entropy - store_r * log(store_r);
			}
		}
		for (int k = 0; k < 28; k++)
		{
			float store_l = float(left_circle[i][k])/float(particles);
			float store_r = float(right_circle[i][k])/float(particles);
			if(store_l > 0)
			{
				entropy = entropy - store_l * log(store_l);
			}
			if(store_r > 0)
			{
				entropy = entropy - store_r * log(store_r);
			}
		}
	}
	
}

void System::Entropy_print_file(string file)
{
	print_file_vector(file, entropys);
}


// print file!
void System::System_print_file(string file)
{
	print_file(file, particles, points);

}



void System::System_print_file_low(string file)
{
	vector<Point> k;
	int store = 0;
	for (int i = 0; i < particles; i++)
	{
		int spot_x = int(points[i].x/.2);
		int spot_y = int(points[i].y/.1);
		if (spot_x == 1 && spot_y == 1)
		{
			k.push_back(start_points[i]);
			store++;
		}
	}
	print_file(file, store, k);

}


//////////////////////////////////////////////////////////////////////
//                       move particle around
/////////////////////////////////////////////////////////////////////

/*        
    ----------              -|
  /            \             | 1
 /              \            | 
|                |          -|
 \              /
  \            /
    ----------
*/


Point after_time_single( Point p, float t)
{
	// new point
	float x = 0;
	float y = 0;
	float v_x = 0;
	float v_y = 0;


	// four times to compair
	float time_upper = (1-p.y)/p.v_y;
	if(time_upper <= 0.00001 || p.v_y == 0)
	{
		time_upper = 10000;
	}

	float time_lower = (-1-p.y)/p.v_y;
	if(time_lower <= 0.00001 || p.v_y == 0)
	{
		time_lower = 10000;
	}
	// left guy
	float a_l = (p.v_x*p.v_x) + (p.v_y*p.v_y);
		
	float b_l = (p.x*p.v_x) + (p.y*p.v_y) + p.v_x;
	
	float c_l = sqrt(-(p.v_x*p.v_x * p.y*p.y) - (p.v_y*p.v_y * p.x*p.x) - (2 * p.v_y*p.v_y * p.x) + (2*(p.v_x*p.v_y*p.x + p.v_x*p.v_y)*p.y) + (p.v_x*p.v_x));

	float time_left_1 = -(b_l-c_l)/a_l;
	
	float time_left_2 = -(b_l+c_l)/a_l;
	
	///printf("a_l = %g, b_l = %g, c_l = %g", a_l, b_l, c_l); 
	
	// error checking
	//print_point(p);
	//printf("time %g\n", t);
	
	
	
	
	
	
	
	
	// new left points
	float x_l_1 = p.x + time_left_1*p.v_x;
	float x_l_2 = p.x + time_left_2*p.v_x;
	// time to use
	float time_left = 10000;
	// check to see if it hits the circle
	if (c_l == 0)
	{
		time_left = 10000;
	}
	// check to see each intersection
	if (time_left_1 > 0.00001 && x_l_1 < -1)
	{
		time_left = time_left_1;
	}
	if (time_left_2 > 0.00001 && x_l_2 < -1)
	{
		time_left = time_left_2;
	}
	
	// right guy
	float b_r = (p.x*p.v_x) + (p.y*p.v_y) - p.v_x;
	
	float c_r = sqrt(-(p.v_x*p.v_x * p.y*p.y) - (p.v_y*p.v_y * p.x*p.x) + (2 * p.v_y*p.v_y * p.x) + (2*(p.v_x*p.v_y*p.x - p.v_x*p.v_y)*p.y) + (p.v_x*p.v_x));
	
	float time_right_1 = -(b_r-c_r)/a_l;
	
	float time_right_2 = -(b_r+c_r)/a_l;
		
	// error checking
	//printf("time right %g, %g \n", time_right_1, time_right_2);
	
	// new right points
	float x_r_1 = p.x + time_right_1*p.v_x;
	float x_r_2 = p.x + time_right_2*p.v_y;
	// time to use
	float time_right = 10000;
	// check to see if it hits the circle
	if (c_r == 0)
	{
		time_left = 10000;
	}
	// check to see each intersection
	if (time_right_1 > 0.00001 && x_r_1 > 1)
	{
		time_right = time_right_1;
	}
	if (time_right_2 > 0.00001 && x_r_2 > 1)
	{
		time_right = time_right_2;
	}

	//printf("time_left = %g, time_right = %g, time_upper = %g, time_lower = %g \n", time_left, time_right, time_upper, time_lower);
	// which side where 1 is top 2 is bottom 3 is left 4 is right 
	int which_side = 0;

	


	if ( t <= time_lower  && t <= time_upper  && t <= time_left  && t <= time_right )
	{
		Point R = {p.x + t*p.v_x, p.y + t*p.v_y, p.v_x, p.v_y};
		return R;
	}


	if (time_upper >= 0.00001 && time_upper <= time_lower && time_upper <= time_left && time_upper <= time_right )
	{
		which_side =  1;
		Point R = {p.x + p.v_x*time_upper, p.y + p.v_y*time_upper, p.v_x, -p.v_y};
		return after_time_single(R, t-time_upper);
	}
	if (time_lower >= 0.00001 && time_lower <= time_upper  && time_lower <= time_left  && time_lower <= time_right )
	{
		which_side =  2;
		Point R = {p.x + p.v_x*time_lower, p.y + p.v_y*time_lower, p.v_x, -p.v_y};
		return after_time_single(R, t-time_lower);
	}
	if ((time_left >= 0.00001 && time_left <= time_upper  && time_left <= time_lower  && time_left <= time_right))
	{
		float c_x = p.x + p.v_x*time_left;
		float c_y = p.y + p.v_y*time_left;
		float n_x = -(c_x+1);
		float n_y = -c_y;
		float u_x = (p.v_x*n_x + p.v_y*n_y)*n_x;
		float u_y = (p.v_x*n_x + p.v_y*n_y)*n_y;
		float w_x = p.v_x - u_x;
		float w_y = p.v_y - u_y;
		Point R = {c_x, c_y, w_x - u_x, w_y - u_y};
		return after_time_single(R, t - time_left);
	}
	if ((time_right >= 0.00001 && time_right <= time_lower  && time_right <= time_left  && time_right <= time_upper))
	{
		float c_x = p.x + p.v_x*time_right;
		float c_y = p.y + p.v_y*time_right;
		float n_x = -(c_x-1);
		float n_y = -c_y;
		float u_x = (p.v_x*n_x + p.v_y*n_y)*n_x;
		float u_y = (p.v_x*n_x + p.v_y*n_y)*n_y;
		float w_x = p.v_x - u_x;
		float w_y = p.v_y - u_y;
		Point R = {c_x, c_y, w_x - u_x, w_y - u_y};
		return after_time_single(R, t - time_right);
	}
	// upper top of the table





}







// ----------------------------------------------------------------------
// Helpful  functions --------------------------------------------------
// ----------------------------------------------------------------------

void print_point(Point p)
{
	printf("(%g, %g) velocity %g, %g \n", p.x, p.y, p.v_x, p.v_y);
}

void print_vector_vector(vector<vector<int> > x)
{
	for(int i = 0; i < x.size(); i++)
	{
		for(int j = 0; j < x[i].size(); j++)
		{
			printf("%d ", x[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}


///////////////////////////////////////////////////////////////
//                    file maker                           ////
//////////////////////////////////////////////////////////////

// makes random files all in lowest box thingy
void print_file_rand_low(string file, int particles)
{
	// seed dis shit
	srand(seed);

	// Make da stream
	ofstream writing;
	writing.open(file.c_str());

	// Write it all down
	writing << particles << "\n";

	for (int i = 0; i < particles; i++)
	{
		writing << (float)rand()/((float)RAND_MAX*50) << " " << (float)rand()/((float)RAND_MAX * 50) << "\n";
		writing << (float)rand()/((float)RAND_MAX*5)-.1 << " " << (float)rand()/((float)RAND_MAX*5)-.1 << "\n";
	}

	// close dis shit
	writing.close();



}


// prints a vector to a file
void print_file_vector(string file, vector<float> k)
{
	// Make da stream
	ofstream writing;
	writing.open(file.c_str());

	// Write it all down

	for (int i = 0; i < k.size(); i++)
	{
		writing << k[i] << "\n";
	}

	// close dis shit
	writing.close();




}


// prints this to file (used in the System stuff)
void print_file(string file, int particles, vector<Point> k)
{

	// Make da stream
	ofstream writing;
	writing.open(file.c_str());

	// Write it all down
	writing << particles << "\n";

	for (int i = 0; i < particles; i++)
	{
		writing << k[i].x << " " << k[i].y << "\n";
		writing << k[i].v_x << " " << k[i].v_y << "\n";
	}

	// close dis shit
	writing.close();

}

// makes random file in all the boxes
void print_file_rand_high(string file, int size, int particles)
{
	// seed dis shit
	srand(seed);

	// Make da stream
	ofstream writing;
	writing.open(file.c_str());

	// Write it all down
	writing << size << "\n" << particles << "\n";

	for (int i = 0; i < particles; i++)
	{
		writing << (float)rand()/((float)RAND_MAX) << " " << (float)rand()/((float)RAND_MAX) << "\n";
	}

	// close dis shit
	writing.close();

}

// prints a universal file, nice for to call in other funcitons


void print_file_rand(string file, int particles)
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




