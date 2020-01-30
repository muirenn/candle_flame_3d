#pragma once

#include "VolumeInstance.h"
#include "Math_Solver.h"

#define DEFAULT_time 30

class FullAnimation
{
private:
	VolumeInstance frame;
	VolumeInstance frame_prev;
	
	double* 
		qty_to_display; // temperature?
	
	unsigned int 
		time_, 
		lin_itr_;
	
	unsigned long int 
		n_;
	
	unsigned long long int 
		frame_size_, 
		iter; //max: ; cubic root [n_] = 2 642 246
	
	double
		dens_f_,
		dens_h_, // density of fuel and hot gaseous products
		Vf_,
		Vh_, // normal velocity of fuel and hot gaseous products
		minTotal_,
		maxTotal_,
		S = 0.161, // reaction speed, g/min; in source - 0.5 m/s; maybe 27.5mm/30min
		ignition_temp_K_ = 472; // or 643?

	string filename_;
	
	void write_init_to_file() const;

public:
	void init();
	void run();
	void find_new_velocity_field(); 
	void write_animation_to_file() const;
	void add_dens_to_center();
	void apply_gravity();
	void resolve_temperature();
	void write_results();
	FullAnimation();
	FullAnimation(int time, int n, string filename);
	~FullAnimation();
};
