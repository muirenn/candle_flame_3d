#pragma once

#include "VolumeInstance.h"

#define DEFAULT_time 30

class FullAnimation
{
private:
	// VolumeInstance* frames;
	VolumeInstance* frame;
	VolumeInstance* frame_prev;
	
	double* dens_volumes;
	unsigned int time_;
	unsigned long int n_;
	unsigned long long int frame_size_, iter; //max: ; cubic root [n_] = 2 642 246
	double minTotal_, maxTotal_;
	string filename_;
	void write_init_to_file() const;

	void project(double* u, double* v, double* w, double* p, double* div) const;
	void set_bnd(int type, double* values) const;
	void diffuse(double dt, double diff, int type, double* values, double* values_0) const;
	void advect(double dt, int type, double* d, double* d0, double* u, double* v, double* w);
	void lin_solve(int type, double* values, double* values_0, double a, double c) const;
	
public:
	void init();
	void run();
	void write_animation_to_file() const;
	FullAnimation();
	FullAnimation(int time, int n, string filename);
	~FullAnimation();
};
