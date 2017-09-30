FROM gcc:latest

WORKDIR /app

# Copy the current directory contents into the container at /app
ADD . /app