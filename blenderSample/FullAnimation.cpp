#include "stdafx.h"
#include "FullAnimation.h"


void FullAnimation::init()
{
	write_init_to_file();

	frames = new VolumeInstance[time_];
	frames[0].allocate_memory(n_);
	 frames[0].draw_sphere();
	// frames[0].draw_candle_v1();

	minTotal_ = frames[0].get_min();
	maxTotal_ = frames[0].get_max();
}

void FullAnimation::run()
{
	for (int i = 1; i < time_; i++)
	{
		frames[i].allocate_memory(n_);
		frames[i].set_prev_values_nochange(frames[i - 1]);
		frames[i].step();

		if (frames[i].get_min() < minTotal_)
		{
			minTotal_ = frames[i].get_min();
		}
		if (frames[i].get_max() > maxTotal_)
		{
			maxTotal_ = frames[i].get_max();
		}
	}
}

void FullAnimation::write_init_to_file() const
{
	ofstream blender_file;
	blender_file.open(filename_, ios::out | ios::binary | ios::trunc);
	blender_file.write(reinterpret_cast<const char*>(&n_), sizeof(n_));
	blender_file.write(reinterpret_cast<const char*>(&n_), sizeof(n_));
	blender_file.write(reinterpret_cast<const char*>(&n_), sizeof(n_));
	blender_file.write(reinterpret_cast<const char*>(&time_), sizeof(time_));
	blender_file.close();
}

void FullAnimation::write_animation_to_file() const
{
	ofstream blender_file;
	blender_file.open(filename_, ios::out | ios::binary | ios::app);

	for (int l = 0; l < time_; l++)
	{
		for (int i = 1; i <= n_; i++)
		{
			for (int j = 1; j <= n_; j++)
			{
				for (int k = 1; k <= n_; k++)
				{
					// C++ writes files in native format. On a standard x86-derived PC native is little endian
					float currentV = static_cast<float>(frames[l].get_density(IDX(i, j, k, n_)) - minTotal_) / (
						maxTotal_ - minTotal_);
				//	float currentV = static_cast<float>(frames[l].get_density(IDX(i, j, k, n_)));
					blender_file.write(reinterpret_cast<const char*>(&currentV), sizeof(currentV));
				}
			}
		}
	}

	blender_file.close();
}

FullAnimation::FullAnimation(): frames(nullptr), time_(DEFAULT_time), n_(DEFAULT_cube_side_size), minTotal_(0), maxTotal_(0)
{
}

FullAnimation::FullAnimation(const int time, const int n, const string filename) : FullAnimation()
{
	time_ = time;
	n_ = n;
	filename_ = filename;
}


FullAnimation::~FullAnimation()
{
	//for(int i = 0; i < time_; i++)
	//{
	//	frames[i].~VolumeInstance();
	//}
	delete[] frames;
}
