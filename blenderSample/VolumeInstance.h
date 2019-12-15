#pragma once

#include "stdafx.h"
#include <iostream>
#include <fstream>
using namespace std;


#define DEFAULT_cube_side_size 50
#define DEFAULT_viscosity 0.0000001f
#define DEFAULT_diffusivity 0
#define DEFAULT_time_step 0.05
#define DEFAULT_avg_density 5
#define DEFAULT_high_density 20
#define DEFAULT_avg_velocity 0.8
#define DEFAULT_high_velocity 1.2
#define DEFAULT_iter 20

constexpr const char* bvox_filename_DEFAULT = "holyshit3.bvox";

#define IDX(i,j,k,N) ((i)+j*((N)+2)+(k)*((N)+2)*((N)+2))
#define SWAP(x0, x) {double *tmp=x0;x0=x;x=tmp;}

class VolumeInstance
{
private:
	enum property { density = 0, velocity_x = 1, velocity_y = 2, velocity_z = 3 };

	int n_;
	unsigned int size_;
	double visc_, diff_, dt_;
	double min_, max_;

	double* density_;
	double* density_prev_;
	double *velocity_x_, *velocity_y_, *velocity_z_;
	double *velocity_x_prev_, *velocity_y_prev_, *velocity_z_prev_;

	void add_source(double* x, double* s) const;
	void project(double* u, double* v, double* w, double* p, double* div) const;
	void set_bnd(int b, double* x) const;
	void diffuse(int b, double* x, double* x0);
	void advect(int b, double* d, double* d0, double* u, double* v, double* w);
	void vel_step();
	void dens_step();

public:
	VolumeInstance();
	~VolumeInstance();
	void allocate_memory();
	void step();
	void draw_sphere() const;
	double get_density(int i) const;
	double get_min() const;
	double get_max() const;
	void set_prev_values(const VolumeInstance& prev_frame);
};
