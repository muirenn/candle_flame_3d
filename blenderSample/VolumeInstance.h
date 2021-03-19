#pragma once

#include "predefined_constants.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <ctime>
using namespace std;

class VolumeInstance
{
private:

public:

	VolumeInstance();
	~VolumeInstance();
	double // 3D arrays
		*dens, // ρ
		*vX,
		*vY,
		*vZ;

	double
		visc_, // ν
		diff_, // κ
		dt_,
		min_,
		max_;
			
	void allocate_memory(unsigned long long int size);
	void add_source(double* x, double* a, unsigned long int n);
	void add_constant(double* x, double a, unsigned long int n);
	void fill_with_zeros(unsigned long int n);
	void hardcode_init_values(unsigned long int n);
	void copy(unsigned long int n, VolumeInstance* from_volume);
	void add_fuel(unsigned long int n);
};
