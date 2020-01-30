#pragma once

#include "stdafx.h"
#include "predefined_constants.h"
// #include "Voxel.h"
#include <fstream>
#include <iostream>
using namespace std;

// The reason it doesn’t work is that the it can become unstable
// when the diffusion rate is too high, time step too large or
// grid spacing too small. The problem is that this method breaks
// down when the density propagates further than just between neighbors.


#define DEFAULT_cube_side_size 50
#define DEFAULT_viscosity 0.001
#define DEFAULT_diffusivity 0 // L^2*T^(-1); L = initial size of gas blob or field
#define DEFAULT_time_step 0.2 // T
#define DEFAULT_iter 20

#define DEFAULT_avg_density 5
#define DEFAULT_high_density 20
#define DEFAULT_avg_velocity 0.2
#define DEFAULT_high_velocity 0

constexpr const char* bvox_filename_DEFAULT = "holyshit.bvox";

#define IDX(i,j,k,N) ((k)+(j)*((N)+2)+(i)*((N)+2)*((N)+2)) // i === z


class VolumeInstance
{
private:

public:

	VolumeInstance();
	~VolumeInstance();
	double // 3D arrays
		*dens_, // ρ
		*vX_,
		*vY_,
		*vZ_,
		*tempK_,
		*press_, // p
		*f_conf,
		*Y_elaps; //1-Y = time elapsed since a fluid element crossed over the blue reaction core; Y(0)=1 in the region occupied by gaseous fuel

	double *f_buoy; // 1D array - Z axis only


	/*
	Positive in the region of space filled with fuel,
	negative elsewhere and zero at the reaction zone
	*/
	int* impl_surf_; 

	//double* unit_normal_;

	double 
		visc_, // ν
		diff_, // κ
		dt_,
		min_,
		max_;
			
	void allocate_memory(unsigned long long int size);
	void add_source(double* x, double* a, unsigned long int n);
	void add_constant(double* x, double a, unsigned long int n);
	void draw_sphere(unsigned long int n);
	void draw_candle_v1(unsigned long int n) const;
	double get_min() const;
	double get_max() const;
	void fill_with_zeros(unsigned long int n);
	void copy(unsigned long int n, VolumeInstance* from_volume);
	void find_f_buoy(unsigned long int n);
	void add_fuel(unsigned long int n);
};
