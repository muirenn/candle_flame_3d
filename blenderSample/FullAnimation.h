#pragma once

#include "VolumeInstance.h"

class FullAnimation
{
private:
	VolumeInstance* frames;

	int time_, n_;
	double minTotal_, maxTotal_;
	string filename_;
	void write_init_to_file() const;

public:
	void init();
	void run();
	void write_animation_to_file() const;
	FullAnimation();
	FullAnimation(int time, int n, string filename);
	~FullAnimation();
};
