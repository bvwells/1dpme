# 1dpme
[![Build Status](https://travis-ci.org/bvwells/1dpme.svg?branch=master)](https://travis-ci.org/bvwells/1dpme)

One dimension Porous Medium Equation (PME)

## Developing

Developing locally requires Docker for Windows. Run the command

```
docker build -t 1dpme .
```

to build the docker image which pulls the gcc base image containing gfortran and maps the source code into the container.

Then run image with the command:

```
docker run -i -t -v /f/git/src/github.com/bvwells/1dpme:/app 1dpme```

This command maps the local workspace into the running image so any changes made in the running image will be reflected on the local workspace.

Within the running image generate the make files by running the command:

```
cmake .
```

Build the executable by running the command:

```
make
```


