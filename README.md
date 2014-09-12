# Billiard Ball Simulator

This program is designed to run large scale billiard ball simulations. It is capable of handling around 10,000,000 balls. It was used to investigate dynamical systems with entropy boundary conditions and perturbations. There is a strange effect that given time boundary conditions the effect of perturbations will follow the direction of entropy increase. It is hypothesized that causality is a result of entropy increase. More information on such matter can be found in my professors paper Causality is an Effect.  http://arxiv.org/abs/cond-mat/0011507

## Text file example

Here is what a text file containing ball information might look like

```
10000 (number of balls

-1  0 (x, y position 
.4  0.0  (velocity_x, velocity_y) 
...

...
```

When running this file using the graphics it produces this simulation



![ScreenShot](https://github.com/loliverhennigh/Bill/blob/master/run_pic.png)


## How it works

The algorithm is very simple. It just iterates through all the balls calculating there final position one at a time. 

## Commands

Here is a list of commands

```
t = timestep of .001

f = timestep of .01

g = timestep of .1

- = reverse direction

1 = fullscreen

2 = exit fullscreen

arrows = rotate

w,a,s,d = move the image
```


## How to compile


To use this program without graphics it is very simple. A example main file called bill.cpp is provided to show all file io commands.

The graphics side is very simple as well. You will need gsl and glut (if you would like graphics). To install these on ubunto simply type

```
sudo apt-get install libgsl0-dev
```

and


```
sudo apt-get install freeglut3-dev
```

Next compile the program

```
make
```

and to execute

```
./bill
```

If glut gives you any trouble consult http://kiwwito.com/installing-opengl-glut-libraries-in-ubuntu/









