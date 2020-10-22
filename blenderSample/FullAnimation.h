#pragma once

#include "VolumeInstance.h"
#include "Math_Solver.h"


class FullAnimation
{
private:
	VolumeInstance frame;
	VolumeInstance frame_prev;
	
	double* 
		qty_to_display;
	
	unsigned int 
		time_, 
		lin_itr_;
	
	unsigned long int 
		n_;
	
	unsigned long long int 
		frame_size_, 
		iter; //max: ; cubic root [n_] = 2 642 246
	
	double
		minTotal_,
		maxTotal_;

	string filename_;
    string log_filename_;

    void write_init_data_to_file() const;

public:
	void init();
	void run();
	void write_animation_to_file() ;
	void apply_buoyancy_forces();
	void apply_confinement_forces();
	void write_results();
	FullAnimation();
	FullAnimation(unsigned int time, unsigned int n, string filename, string log_filename);
	~FullAnimation();
};
