#include "stdafx.h"
#include "FullAnimation.h"

void FullAnimation::init()
{
	iter = 0;
	write_init_to_file();

	frame_size_ = (n_ + 2) * (n_ + 2) * (n_ + 2);

	frame.allocate_memory(frame_size_);
	frame_prev.allocate_memory(frame_size_);

	unsigned long long int d_length = time_ * frame_size_;
	qty_to_display = new (nothrow) double[d_length];
	if(qty_to_display == nullptr)
	{
		cerr << "MEMORY ALLOCATION ERROR" << endl;
	}
	
	frame.draw_sphere(n_);
	//frame.draw_candle_v1(n_);
	frame_prev.fill_with_zeros(n_);

	minTotal_ = T_MAX;
	maxTotal_ = -20;
}

void FullAnimation::run()
{
	// define/track the blue core


	unsigned int t;

	for (t = 0; t < time_; t++)
	{
		/* Velocity step: ∂u/∂t=-(u∙∇)u+ν∇^2 u+f */
		//apply_buoyancy_forces();

		diffuse(frame.dt_, frame.diff_, velocity_x, frame_prev.vX_, frame.vX_, lin_itr_, n_);
		diffuse(frame.dt_, frame.diff_, velocity_y, frame_prev.vY_, frame.vY_, lin_itr_, n_);
		diffuse(frame.dt_, frame.diff_, velocity_z, frame_prev.vZ_, frame.vZ_, lin_itr_, n_);

		project(frame_prev.vX_, frame_prev.vY_, frame_prev.vZ_, frame.vX_, frame.vY_, lin_itr_, n_);

		advect(frame.dt_, velocity_x, frame.vX_, frame_prev.vX_, frame_prev.vX_, frame_prev.vY_, frame_prev.vZ_, n_);
		advect(frame.dt_, velocity_y, frame.vY_, frame_prev.vY_, frame_prev.vX_, frame_prev.vY_, frame_prev.vZ_, n_);
		advect(frame.dt_, velocity_z, frame.vZ_, frame_prev.vZ_, frame_prev.vX_, frame_prev.vY_, frame_prev.vZ_, n_);

		project(frame.vX_, frame.vY_, frame.vZ_, frame_prev.vX_, frame_prev.vY_, lin_itr_, n_);

		// TODO add pressure term 
		// solve ∇∙(∇p/ρ)=∇∙u^*/∆t


		/* Density step: ∂ρ/∂t=-(u∙∇)ρ+κ∇^2 ρ+S*/

		//diffuse(frame.dt_, frame.diff_, density, frame_prev.dens_, frame.dens_, lin_itr_, n_);
		//advect(frame.dt_, density, frame.dens_, frame_prev.dens_, frame.vX_, frame.vY_, frame.vZ_, n_);


		/* Reaction coordinate variable step */
		diffuse(frame.dt_, frame.diff_, density, frame_prev.Y_elaps, frame.Y_elaps, lin_itr_, n_);

		advect(frame.dt_, density, frame.Y_elaps, frame_prev.Y_elaps, frame.vX_, frame.vY_, frame.vZ_, n_);
		frame.add_constant(frame.Y_elaps, -CONST_K, n_);
		//frame.add_constant(frame.vY_, DEFAULT_avg_velocity, n_);

		frame_prev.copy(n_, &frame);

		/* Temperature step */
		//resolve_temperature();




		write_results();

		iter++;

		//frame_prev.copy(n_, &frame);
		//add_dens_to_center();
		//frame.add_fuel(n_);
	}
}

void FullAnimation::find_new_velocity_field()
{








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

	unsigned long int i, j, k, l;

	
	for(l = 0; l < time_; l++)
	{
		for (k = 1; k <= n_; k++)
		{
			for (j = 1; j <= n_; j++)
			{
				for (i = 1; i <= n_; i++)
				{
					float qty_normal = static_cast<float>((qty_to_display[l*frame_size_ + IDX(i, j, k, n_)] - minTotal_) / (maxTotal_ - minTotal_));
					//float qty_normal = static_cast<float>(qty_to_display[l*frame_size_ + IDX(i, j, k, n_)]);
					blender_file.write(reinterpret_cast<const char*>(&qty_normal), sizeof(qty_normal));
				}
			}
		}
	}

	
	blender_file.close();
}

void FullAnimation::add_dens_to_center()
{
	long int i, j, k;

		long int half_n = static_cast<long int>(n_ / 2);
		for (i = -2; i <= 2; i++)
		{
			for (j = -2; j <= 2; j++)
			{
				for (k = -2; k <= 2; k++)
				{
					int ad = 9 - abs(i*j*k);
					//cout << ad << endl;
					frame.dens_[IDX(half_n + i, half_n + j, half_n + k, n_)] += ad * 100.0;
				}
			}
		}
		double vx = static_cast <double> (rand()) / static_cast <double> (RAND_MAX) - 0.5;
		double vy = static_cast <double> (rand()) / static_cast <double> (RAND_MAX) - 0.5;
		double vz = static_cast <double> (rand()) / static_cast <double> (RAND_MAX) - 0.5;
		//cout << vx << " " << vy << " " << vz << endl;
		frame.vX_[IDX(half_n, half_n, half_n, n_)] += vx * 10;
		frame.vY_[IDX(half_n, half_n, half_n, n_)] += vy * 10;

		frame.vZ_[IDX(half_n, half_n, half_n, n_)] += vz * 10;

	
}

void FullAnimation::apply_buoyancy_forces()
{
		unsigned long int i, j, k;

		for (i = 1; i <= n_; i++)
		{
			for (j = 1; j <= n_; j++)
			{
				for (k = 1; k <= n_; k++)
				{
					//cout << ad << endl;
					frame.vZ_[IDX(i, j, k, n_)] += ALPHA * (frame.tempK_[IDX(i, j, k, n_)] - KELVINIZE(T_AIR));
				}
			}
		}
}

void FullAnimation::resolve_temperature()
{
	unsigned long int i, j, k;

	for (i = 1; i <= n_; i++)
	{
		for (j = 1; j <= n_; j++)
		{
			for (k = 1; k <= n_; k++)
			{
				const double y = frame.Y_elaps[IDX(i, j, k, n_)];
				if (y <= 1 && y >= 0.9) {
					frame.tempK_[IDX(i, j, k, n_)] = KELVINIZE(y * 10 * (T_IGN - T_MAX) + 10 * T_MAX - 9 * T_IGN);
				}
				else {
					if (y > 1) {
						frame.tempK_[IDX(i, j, k, n_)] = KELVINIZE(-y * 10 * (T_IGN - T_MAX) - 10 * T_MAX + 9 * T_IGN);
					//	frame.Y_elaps[IDX(i, j, k, n_)] = 0;
					}

				}
				//cout << UNKELVINIZE(frame.tempK_[IDX(i, j, k, n_)]) << endl;
			}
		}
	}

}

void FullAnimation::write_results()
{
	unsigned long int i, j, k;


	for (i = 1; i <= n_; i++)
	{
		for (j = 1; j <= n_; j++)
		{
			for (k = 1; k <= n_; k++)
			{
				const double res = frame.Y_elaps[IDX(i, j, k, n_)];
				qty_to_display[iter*frame_size_ + IDX(i, j, k, n_)] = res; // TODO check order
				//iter++;
				if (res < minTotal_)
				{
					minTotal_ = res;
				}
				if (res > maxTotal_)
				{
					maxTotal_ = res;
				}	
			}
		}
	}
}

FullAnimation::FullAnimation(): qty_to_display(nullptr), time_(DEFAULT_time),
                                n_(DEFAULT_cube_side_size),
								lin_itr_(LIN_ITERS),
                                frame_size_(0), iter(0),
                                minTotal_(0),
                                maxTotal_(0)
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
	delete[] qty_to_display;
}
