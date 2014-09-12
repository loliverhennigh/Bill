//
// File:        cude.cpp
// Author:      Oliver Hennigh
// Created:     6/1/2013
// Project:     Diffution simulator
// Description: Creates an OpenGL window and draws room
//              That the user can rotate using the arrow keys
// 
// Controls:    Left Arrow  - Rotate Left
//              Right Arrow - Rotate Right
//              Up Arrow    - Rotate Up
//              Down Arrow  - Rotate Down     
//
// ----------------------------------------------------------
// Includes
// ----------------------------------------------------------
//standerd stuff
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <math.h>
#define GL_GLEXT_PROTOTYPES
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

// gsl stuff
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_complex.h>	    
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

// open gl stuff
#pragma comment (lib,"opengl32.lib")
#pragma comment (lib,"glu32.lib")
#pragma comment (lib,"glut32.lib")      

// my stuff
#include "bill.h"

using namespace std;



double scale = 10;
 System test("test.txt");
// ----------------------------------------------------------
// Function Prototypes
// ----------------------------------------------------------
void display();
void specialKeys();
string double_to_string(double number);
 
// ----------------------------------------------------------
// Global Variables
// ----------------------------------------------------------
double rotate_y=0; 
double rotate_x=0;
double xtranslate = 0;
double ytranslate = 0;
double ztranslate = 0;

int width = 200;
int height = 400;

int lastX = 150;
int lastY = 150;

double scale_change = .5;

bool in_full_screen = 0;

// time steps -----------------------------------------------
double dt = .001;
double time_elapsed = 0.0;
int direction = 1;
//  ---------------------------------------------------------

string file_name;

void drawStrokeText(char*string,int x,int y,int z);

void drawText(const char * message, double x, double y, double z);

string heat = "temp ";
string timey = "time = 0";
const char *mess = timey.c_str();
const char *temp = heat.c_str();
// ----------------------------------------------------------
// display() Callback function
// ----------------------------------------------------------
void display(){
 //  Clear screen and Z-buffer
 
  // Reset transformations
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  glColor3f(1,1,1);
  drawText(mess,-.9,-.9,-1);
  glLoadIdentity();    
  // Rotate when user changes rotate_x and rotate_y
  glRotatef( rotate_x, 1.0, 0.0, 0.0 );
  glRotatef( rotate_y, 0.0, 1.0, 0.0 );
  // Movement
  glTranslatef(xtranslate,ytranslate,ztranslate);
  // Scale stuff
  glScalef(scale_change,scale_change,scale_change);
  
  glPointSize(3.0);




  //Multi-colored side - FRONT
  glPushMatrix();

  glPushMatrix();
	glBegin(GL_POLYGON);
        	glColor3f(1, 0, 0); glVertex3f(-1, -1, 0);
        	glColor3f(1, 0, 0); glVertex3f(1, -1, 0);
		glColor3f(1, 0, 0); glVertex3f(1, 1, 0);
		glColor3f(1, 0, 0); glVertex3f(-1, 1, 0);
	glEnd();
  
  glPopMatrix();

  glPushMatrix();
	glBegin(GL_TRIANGLE_FAN);
		glColor3f(0.0f, 0.0f, 1.0f);  // Blue
        	glVertex3f(-1.0, 0.0f, .1);       // Center of circle
	        int numSegments = 100;
		GLfloat angle;
		for (int i = 0; i <= numSegments; i++) { 
			// Last vertex same as first vertex
	                angle = i * 2.0f * 3.1415 / numSegments;  
			// 360 deg for all segments
			glVertex3f(cos(angle)-1, sin(angle), .1);
		}
	glEnd();
  glPopMatrix();
	
  glPushMatrix();
	glBegin(GL_TRIANGLE_FAN);
		glColor3f(0.0f, 0.0f, 1.0f);  // Blue
        	glVertex3f(1.0, 0.0f, .1);       // Center of circle
		for (int i = 0; i <= numSegments; i++) { 
			// Last vertex same as first vertex
	                angle = i * 2.0f * 3.1415 / numSegments;  
			// 360 deg for all segments
			glVertex3f(cos(angle)+1, sin(angle), .1);
		}
	glEnd();
  glPopMatrix();
	
 // draw the nodes (and heats for them)
   for (int i = 0; i < test.particles; i++)
   {
  	glPushMatrix();
		glTranslatef(test.points[i].x,test.points[i].y,0);
		glScalef(.05,.05,.05);
		glColor3f(1,1,1);
  		glutWireTeapot(1);
  	glPopMatrix();
   }

  glFlush();
  glutSwapBuffers();
 
}


void drawText(const char * message, double x, double y, double z)
{
	/* raster pos sets the current raster position
	 * mapped via the modelview and projection matrices
	 */
	glRasterPos3d((GLdouble)x, (GLdouble)y, (GLdouble)z);
	glColor3f(1,1,1);

	// write using bitmap

	while (*message) {
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *message++);
	}
}

 


// ----------------------------------------------------------
// specialKeys() Callback Function
// ----------------------------------------------------------
void specialKeys( int key, int x, int y ) {
 
  //  Right arrow - increase rotation by 5 degree
  if (key == GLUT_KEY_RIGHT)
    rotate_y += 5;
 
  //  Left arrow - decrease rotation by 5 degree
  else if (key == GLUT_KEY_LEFT)
    rotate_y -= 5;
 
  else if (key == GLUT_KEY_UP)
    rotate_x += 5;
 
  else if (key == GLUT_KEY_DOWN)
    rotate_x -= 5;
 
  //  Request display update
  glutPostRedisplay();
 
}

void timeMaker(unsigned char key, int x, int y)
{
	// change time stuff
	if (key ==  't')
	{
		time_elapsed =+ dt*direction;
		test.after_time(dt*direction);
		timey = "time = ";
		timey.append(double_to_string(time_elapsed));
		
		mess = timey.c_str();
	}
	if (key ==  'p')
	{
		test.preturb_step(1);
		timey = "time = ";
		timey.append(double_to_string(time_elapsed));
		
		mess = timey.c_str();
	}
	if (key == 'f')
	{
		time_elapsed =+ dt*10.0*direction;
		test.after_time(dt*10.0*direction);
		
		timey = "time = ";
		timey.append(double_to_string(time_elapsed));
		mess = timey.c_str();

	}
	if (key == 'g')
	{
		time_elapsed += dt*100.0*direction;
		test.after_time(dt*100.0*direction);
		
		timey = "time = ";
		timey.append(double_to_string(time_elapsed));
		mess = timey.c_str();

	}

	if (key == '-')
	{
		direction = direction*(-1);
	}
	
	// move around

	if(key == 'w') 
	{
		ztranslate += 0.02*cos(rotate_y*M_PI/180);
		xtranslate -= 0.02*sin(rotate_y*M_PI/180);	
	}
	if(key == 's')
	{
		ztranslate -= 0.02*cos(rotate_y*M_PI/180);
		xtranslate += 0.02*sin(rotate_y*M_PI/180);
	}
	if(key == 'a')
	{
		xtranslate += 0.02*cos(rotate_y*M_PI/180);
		ztranslate += 0.02*sin(rotate_y*M_PI/180);	
	}
	if(key == 'd')
	{
		xtranslate -= 0.02*cos(rotate_y*M_PI/180);
		ztranslate -= 0.02*sin(rotate_y*M_PI/180);
	}

	// rotate!

	if (rotate_y > 360) rotate_y -= 360;
	else if (rotate_y < -360) rotate_y += 360;

	if(key == 'e')
	{
		scale_change += .001;
	}
	if(key == 'r')
	{
		scale_change -= .001;
	}

	glutPostRedisplay();

	// FULL SCREEN!

	if(key == '1')
	{
		glutFullScreen();
		in_full_screen = 1;
	}
	if(key == '2')
	{
		glutPositionWindow(0,0);
		in_full_screen = 0;
	}
	
	// Print dis (not working yet)	

	if(key == 'p')
	{
		string p = file_name; 
		p.append("_");
		string str = double_to_string(time_elapsed);
		p.append(str);
		p.append(".txt");
		timey.append("   printed to ");
		timey.append(p);
		mess = timey.c_str();
	}

  glFlush();

}


string double_to_string(double number)
{
	ostringstream s;
	s << number;
	return s.str(); 
}


void prompt()
{
	cout << "Enter the text file you want to read data from \n";
	cin >> file_name;
	cout << "Hope that works, I have no error checking :) \n";
	cout << "if its blank then you probably didnt enter the correct file. Have a Great Day ;) \n";
	System sample(file_name);
	test = sample;
}
 


void drawStrokeText(char*string,int x,int y,int z)
{
          char *c; 
          glPushMatrix();
          glTranslatef(x, y,z);
          glScalef(0.09f,-0.08f,z);
  
          for (c=string; *c != '\0'; c++)
          {   
                glutStrokeCharacter(GLUT_STROKE_ROMAN , *c);
          }   
          glPopMatrix();
}

// ----------------------------------------------------------
// main() function
// ----------------------------------------------------------
int main(int argc, char* argv[]){

  //  Prompt function
  prompt(); 

  //  Initialize GLUT and process user parameters
  glutInit(&argc,argv);
 
  //  Request double buffered true color window with Z-buffer
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
 
  // Create window
  glutCreateWindow("q walk");

  //  Enable Z-buffer depth test
  glEnable(GL_DEPTH_TEST);
 
  glutSetCursor(GLUT_CURSOR_CROSSHAIR);
  // Callback functions
  glutDisplayFunc(display);
  glutSpecialFunc(specialKeys);
  glutKeyboardFunc(timeMaker);
//  glutMotionFunc(mouse_movements); 

  //  Pass control to GLUT for events
  glutMainLoop();
 
  //  Return to OS
  return 0;
 
}
