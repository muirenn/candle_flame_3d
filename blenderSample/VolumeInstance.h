#pragma once

#include "stdafx.h"
#include <iostream>
#include <fstream>
using namespace std;


constexpr auto cube_side_size = 30;
constexpr auto viscosity = 0.0000001f;
constexpr auto diffusivity = 0;
constexpr auto time_step = 0.1;
constexpr auto bvox_filename = "holyshit3.bvox";;
constexpr auto avg_density = 25;
constexpr auto high_density = 100;
constexpr auto avg_velocity = 0.8;
constexpr auto high_velocity = 1.2;
constexpr auto iter = 64;
constexpr auto scale = 5;

constexpr void SWAP(float* x0, float* x)
{
	float* tmp = x0;
	x0 = x;
	x = tmp;
}


class VolumeInstance
{
private:
	enum property { density = 0, velocity_x = 1, velocity_y = 2, velocity_z = 3 };

	int n_;
	unsigned int size_;
	float visc_, diff_, dt_;
	float min_, max_;

	float* density_prev_;
	float *velocity_x_, *velocity_y_, *velocity_z_;
	float *velocity_x_prev_, *velocity_y_prev_, *velocity_z_prev_;


	unsigned int idx(int i, int j, int k) const;
	void add_source(float* x, float* s);
	void project(float* u, float* v, float* w, float* p, float* div);
	void set_bnd(int b, float* x);
	void diffuse(int b, float* x, float* x0);
	void advect(int b, float* d, float* d0, float* u, float* v, float* w);
	void vel_step(float* velocity_x, float* velocity_y, float* velocity_z, float* velocity_x_prev,
	              float* velocity_y_prev, float* velocity_z_prev);
	void dens_step(float* density, float* density_prev, float* velocity_x, float* velocity_y, float* velocity_z);
	void add_prev();
	void find_minmax(float* arr);

public:
	VolumeInstance();
	~VolumeInstance();
	void init();
	void step();
	string filename_;
	float* density_;

	void write_to_file(float* arr);
};
