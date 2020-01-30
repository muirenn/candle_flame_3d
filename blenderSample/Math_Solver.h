#pragma once

#include "stdafx.h"
#include <cmath>

#define IDX(i,j,k,N) ((k)+(j)*((N)+2)+(i)*((N)+2)*((N)+2))

enum property { density, velocity_x, velocity_y, velocity_z };

void lin_solve(property prop, double* x, double* x0, double a, double c, int lin_itr, unsigned long int n);
void set_bnd(property prop, double* x, unsigned long int n);
void project(double* u, double* v, double* w, double* p, double* div, int lin_itr, unsigned long int n);
void advect(double dt, property prop, double* d, double* d0, double* u, double* v, double* w, int n);
void diffuse(double dt, double diff, property prop, double* x, double* x0, int lin_itr, unsigned long int n);
