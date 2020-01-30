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
		*imp_surf_vX,
		*imp_surf_vY,
		*imp_surf_vZ,
		*unit_norm, //local unit normal
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
	double get_min() const;
	double get_max() const;
	void fill_with_zeros(unsigned long int n);
	void copy(unsigned long int n, VolumeInstance* from_volume);
	void add_fuel(unsigned long int n);
};
