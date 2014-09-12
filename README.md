# Billiard Ball Simulator

This program is designed to run large scale billiard ball simulations. It is capable of handling around 10,000,000 balls in



## Text file example

Here is what a text file containing ball information might look like

```
10000 (number of balls

-1  0  .4  0.0  (four numbers, x, y, velocity_x, velocity_y) 

...
```

When running this file it produces this screen


Here is another example of a file called pres.txt

![ScreenShot](https://github.com/loliverhennigh/Quantum-Walk-Simulator/blob/master/pres.png)


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


To use this program without graphics it is very simple. A example main file is provied to show all file io commands.

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

and to exacute

```
./bill
```

If glut gives you any trouble consult http://kiwwito.com/installing-opengl-glut-libraries-in-ubuntu/









