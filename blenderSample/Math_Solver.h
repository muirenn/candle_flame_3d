#pragma once

#include <cmath>
#include <string>
#include <bits/stdc++.h>
#include <filesystem>
#include "predefined_constants.h"
using namespace std;

namespace FS = filesystem;

enum property
{
    density,
    velocity_x,
    velocity_y,
    velocity_z,
    dens2
};

void lin_solve(property prop, double *x, double *x0, double a, double c, int lin_itr, unsigned long int n);
void set_bnd(property prop, double *x, unsigned long int n);
void project(double *u, double *v, double *w, double *p, double *div, int lin_itr, unsigned long int n);
void advect(double dt, property prop, double *d, double *d0, double *u, double *v, double *w, int n);
void diffuse(double dt, double b, property prop, double *x, double *x0, int lin_itr, unsigned long int n);
double curl(double *u, double *v, double *w, unsigned long int i, unsigned long int j, unsigned long int k, unsigned long int n);
bool fuel_zone(long int i, long int j, long int k);

void get_filenames(FS::path *log_filename, FS::path *bvox_filename);