# 1dpme
One dimension Porous Medium Equation (PME)

## Developing locally.

Developing locally requires Docker for Windows. Run the command

```
docker build -t 1dpme .
``

to build the docker image which pulls the gcc base image containing gfortran and maps the source code into the container.

Then run image with the command:

```
docker run -i -t 1dpme
```

