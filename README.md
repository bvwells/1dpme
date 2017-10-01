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
u(x,t) = (1/λ)*(1-(x/(r~0~*λ)))^(1/m)
```

## Developing

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

[1]: http://www.reading.ac.uk/nmsruntime/saveasdialog.aspx?lID=24080&sID=90294
