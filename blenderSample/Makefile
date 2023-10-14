# Makefile for Writing Make Files Example

# *****************************************************
# Variables to control Makefile operation

CC = g++
CFLAGS = -Wall -g

# ****************************************************
# Targets needed to bring the executable up to date

blenderSample: blenderSample.o FullAnimation.o VolumeInstance.o Math_Solver.o
	$(CC) $(CFLAGS) -o blenderSample blenderSample.o FullAnimation.o VolumeInstance.o Math_Solver.o

# The main.o target can be written more simply

blenderSample.o: blenderSample.cpp FullAnimation.h VolumeInstance.h Math_Solver.h
	$(CC) $(CFLAGS) -c blenderSample.cpp

FullAnimation.o: FullAnimation.h VolumeInstance.h Math_Solver.h

VolumeInstance.o: VolumeInstance.h Math_Solver.h predefined_constants.h

Math_Solver.o: Math_Solver.h predefined_constants.h