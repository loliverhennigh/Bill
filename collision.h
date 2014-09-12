#include <math.h>
#include "Points.h"

using namespace std;


Point after_time_single( Point p, float t)
{
	// new point
	float x = 0;
	float y = 0;
	float v_x = 0;
	float v_y = 0;


	// four times to compair
	float time_upper = (1-p.y)/p.v_y;
	if(time_upper <= 0 || p.v_y == 0)
	{
		time_upper = -10000;
	}

	float time_lower = (-1-p.y)/p.v_y;
	if(time_lower <= 0 || p.v_y == 0)
	{
		time_lower = -10000;
	}
	// left guy
	float a_l = p.v_x*p.v_x+p.v_y*p.v_y;
	float b_l = p.x*p.v_x + p.y*p.v_y + p.v_x;
	float c_l = sqrt(- p.v_x*p.v_x * p.y*p.y - p.v_y*p.v_y * p.x*p.x - 2 * p.v_y*p.v_y * p.x + 2*(p.v_x*p.v_y*p.x + p.v_x*p.v_y)*p.y + p.v_x*p.v_x);
	float time_left_1 = -(b_l-c_l)/a_l;
	float time_left_2 = -(b_l+c_l)/a_l;
	// new left points
	float x_l_1 = p.x + time_left_1*p.v_x;
	float x_l_2 = p.x + time_left_2*p.v_x;
	// time to use
	float time_left = -10000;
	// check to see if it hits the circle
	if (c_l == 0)
	{
		time_left = -10000;
	}
	// check to see each intersection
	if (time_left_1 > 0 && x_l_1 < -1)
	{
		time_left = time_left_1;
	}
	if (time_left_2 > 0 && x_l_2 < -1)
	{
		time_left = time_left_2;
	}
	
	// right guy
	float c_r = sqrt(- p.v_x*p.v_x * p.y*p.y - p.v_y*p.v_y * p.x*p.x - 2 * p.v_y*p.v_y * p.x + 2*(p.v_x*p.v_y*p.x + p.v_x*p.v_y)*p.y - p.v_x*p.v_x - 2*p.v_y*p.v_y);
	float time_right_1 = -(b_l-c_r)/a_l;
	float time_right_2 = -(b_l+c_r)/a_l;
	// new right points
	float x_r_1 = p.x + time_right_1*p.v_x;
	float x_r_2 = p.x + time_right_2*p.v_y;
	// time to use
	float time_right = -10000;
	// check to see if it hits the circle
	if (c_r == 0)
	{
		time_left = -10000;
	}
	// check to see each intersection
	if (time_right_1 > 0 && x_r_1 > 1)
	{
		time_right = time_right_1;
	}
	if (time_right_2 > 0 && x_r_2 > 1)
	{
		time_right = time_right_2;
	}

	printf("%g, %g, %g, %g,", time_lower, time_upper, time_left, time_right);
	// which side where 1 is top 2 is bottom 3 is left 4 is right 
	int which_side = 0;
	if ( t < abs(time_lower) && t < abs(time_upper) && t < abs(time_left) && t < abs(time_right))
	{
		Point R = {p.x + t*p.v_x, p.y + t*p.v_y, p.v_x, p.v_y};
		print_point(R);
		return R;
	}


	if (time_upper > 0 && time_upper < abs(time_lower) && time_upper < abs(time_left) && time_upper < abs(time_right))
	{
		which_side =  1;
		Point R = {p.x + p.v_x*time_upper, p.y + p.v_y*time_upper, p.v_x, -p.v_y};
		return after_time_single(R, t-time_upper);
	}
	if (time_lower > 0 && time_lower < abs(time_upper) && time_lower < abs(time_left) && time_lower < abs(time_right))
	{
		which_side =  2;
		Point R = {p.x + p.v_x*time_lower, p.y + p.v_y*time_lower, p.v_x, -p.v_y};
		return after_time_single(R, t-time_lower);
	}
	if (time_left > 0 && time_left < abs(time_upper) && time_left < abs(time_lower) && time_left < abs(time_right))
	{
		which_side =  3;
	}
	if (time_upper > 0 && time_upper < abs(time_lower) && time_upper < abs(time_left) && time_upper < abs(time_right))
	{
		which_side =  4;
	}
	// upper top of the tabel
	





}



