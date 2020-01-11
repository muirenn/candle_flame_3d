#pragma once

#include "stdafx.h"
#include "startup_values.h"
#include <fstream>
#include <iostream>
using namespace std;

// The reason it doesn’t work is that the it can become unstable
// when the diffusion rate is too high, time step too large or
// grid spacing too small. The problem is that this method breaks
// down when the density propagates further than just between neighbors.


#define DEFAULT_cube_side_size 50
#define DEFAULT_viscosity 0.0000001
#define DEFAULT_diffusivity 0 // L^2*T^(-1); L = initial size of gas blob or field
#define DEFAULT_time_step 0.05 // T
#define DEFAULT_iter 4

#define DEFAULT_avg_density 5
#define DEFAULT_high_density 20
#define DEFAULT_avg_velocity 0.5
#define DEFAULT_high_velocity 1.2

constexpr const char* bvox_filename_DEFAULT = "holyshit.bvox";

#define IDX(i,j,k,N) ((i)+j*((N)+2)+(k)*((N)+2)*((N)+2))
#define SWAP(x0, x) {double *tmp=x0;x0=x;x=tmp;}

enum property { density = 0, velocity_x = 1, velocity_y = 2, velocity_z = 3 };

class VolumeInstance
{
private:

public:

	VolumeInstance();
	~VolumeInstance();
	double* density_;
	double *velocity_x_, *velocity_y_, *velocity_z_;
	double visc_, diff_, dt_;
	double min_, max_;
	
	void allocate_memory(unsigned long long int size);
	void step();
	void draw_sphere(unsigned long int n);
	void draw_candle_v1(unsigned long int n) const;
	double get_min() const;
	double get_max() const;
	void fill_with_zeros(unsigned long int n);
};
