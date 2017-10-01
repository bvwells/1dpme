# 1dpme
[![Build Status](https://travis-ci.org/bvwells/1dpme.svg?branch=master)](https://travis-ci.org/bvwells/1dpme)

One dimensional Porous Medium Equation (PME) solved by a moving mesh approach
described in the PhD thesis,

*A moving mesh finite element method for the numerical solution of partial differential equations and systems.*

which can be found [here][1].

The Porous Medium Equation is described by the non-linear partial differential equation

```
u_t = (u^m u_x)_x
```

and admits self-similar solutions of the form

```
u(x,t) = (1/lambda)*(1-(x/(r0*lambda))^2)^(1/m)
```

where

```
r0 = Q*gamma(1/m + 3/2)/(SQRT(pi)*gamma(1/m + 1))                   
t0 = ((r0^2)*m)/(2*(m+2))                            
lambda = (t/t0)^(1/(m+2)) 
```

Here ```Q``` is the total mass of the solution ```gamma``` is the gamma function.

## Numerical Solution

The one-dimensional Porous Medium Equation is solved using a moving mesh 
method which uses the monitor function ```M=u(x,t)``` in the moving mesh 
equations for the mesh velocity. The mesh is advanced forwards in time 
using a forward Euler time-stepping scheme. The solution to the PME 
is obtained as a reconstruction step by considering conservation of
the monitor function. All the moving mesh equations are solved using
linear finite elements.

## Building and Developing

Developing locally requires Docker for Windows. Run the command

```
docker build -t 1dpme .
```

to build the docker image which pulls the gcc base image containing gfortran and maps the source code into the container.

Then run image with the command:

```
docker run -i -t -v /f/git/src/github.com/bvwells/1dpme:/app 1dpme
```

This command maps the local workspace into the running image so any changes made in the running image will be reflected on the local workspace.

Within the running image generate the make files by running the command:

```
cmake .
```

Build the executable by running the command:

```
make
```

## Running

The program takes the file [variables.data](./variables.data) as input to the simulation. The program can be run from the base of the repo with the command:

```
./bin/1dpme.exe
```

The program outputs the mesh and solution over time into the file ```solution.out```.

[1]: http://www.reading.ac.uk/nmsruntime/saveasdialog.aspx?lID=24080&sID=90294
