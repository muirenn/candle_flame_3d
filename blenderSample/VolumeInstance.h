#pragma once

#include "stdafx.h"
#include <iostream>
#include <fstream>
using namespace std;


constexpr auto cube_side_size = 50;
constexpr auto viscosity = 0.0000001f;
constexpr auto diffusivity = 0;
constexpr auto time_step = 0.01;
constexpr auto bvox_filename = "holyshit3.bvox";;
constexpr auto avg_density = 5;
constexpr auto high_density = 20;
constexpr auto avg_velocity = 0.8;
constexpr auto high_velocity = 1.2;
constexpr auto iter = 20;
constexpr auto scale = 5;

//constexpr int idx(int i, int j, int k, int N) { return i+j*(N+2)+k*(N+2)*(N+2); }
#define idx(i,j,k,N)  ((i)+j*((N)+2)+(k)*((N)+2)*((N)+2))


class VolumeInstance
{
private:
	enum property { density = 0, velocity_x = 1, velocity_y = 2, velocity_z = 3 };

	int n_;
	unsigned int size_;
	double visc_, diff_, dt_;
	double min_, max_;

	double* density_prev_;
	double *velocity_x_, *velocity_y_, *velocity_z_;
	double *velocity_x_prev_, *velocity_y_prev_, *velocity_z_prev_;


	// unsigned int idx(int i, int j, int k) const;
	void add_source(double* x, double* s);
	void project(double* u, double* v, double* w, double* p, double* div);
	void set_bnd(int b, double* x);
	void diffuse(int b, double* x, double* x0);
	void advect(int b, double* d, double* d0, double* u, double* v, double* w);
	void vel_step(double* velocity_x, double* velocity_y, double* velocity_z, double* velocity_x_prev,
	              double* velocity_y_prev, double* velocity_z_prev);
	void dens_step(double* density, double* density_prev, double* velocity_x, double* velocity_y, double* velocity_z);
	void add_prev();
	void find_minmax(double* arr);

public:
	VolumeInstance();
	~VolumeInstance();
	void init();
	void step();
	string filename_;
	double* density_;

	void write_to_file(double* arr);
};
