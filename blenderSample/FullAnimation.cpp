#include "stdafx.h"
#include "FullAnimation.h"


void FullAnimation::lin_solve(int type, double* values, double* values_0, double a, double c) const
{
	unsigned long int i, j, k, l;
	for (l = 0; l < DEFAULT_iter; l++)
	{
		for (i = 1; i <= n_; i++)
		{
			for (j = 1; j <= n_; j++)
			{
				for (k = 1; k <= n_; k++)
				{
					values[IDX(i, j, k, n_)] = (values_0[IDX(i, j, k, n_)] +
						a * (
							values[IDX(i - 1, j, k, n_)] +
							values[IDX(i + 1, j, k, n_)] +
							values[IDX(i, j - 1, k, n_)] +
							values[IDX(i, j + 1, k, n_)] +
							values[IDX(i, j, k - 1, n_)] +
							values[IDX(i, j, k + 1, n_)]
							)
						) / c;
				}
			}
		}
		set_bnd(type, values);
	}
}

void FullAnimation::init()
{
	iter = 0;
	write_init_to_file();

	frame = new VolumeInstance();
	frame_prev = new VolumeInstance();

	frame_size_ = (n_ + 2) * (n_ + 2) * (n_ + 2);

	frame->allocate_memory(frame_size_);
	frame_prev->allocate_memory(frame_size_);

	unsigned long long int d_length = time_ * frame_size_;
	dens_volumes = new (nothrow) double[d_length];
	if(dens_volumes == nullptr)
	{
		cout << "MEMORY ALLOCATION ERROR " << endl;
	}
	
	frame->draw_sphere(n_);
	frame_prev->fill_with_zeros(n_);

	minTotal_ = frame->get_min();
	maxTotal_ = frame->get_max();
}

void FullAnimation::run()
{
	for (unsigned int i = 0; i < time_; i++)
	{
		// velocity step
		diffuse(frame->dt_, frame->diff_, velocity_x, frame_prev->velocity_x_, frame->velocity_x_);
		diffuse(frame->dt_, frame->diff_, velocity_y, frame_prev->velocity_y_, frame->velocity_y_);
		diffuse(frame->dt_, frame->diff_, velocity_z, frame_prev->velocity_z_, frame->velocity_z_);

		project(frame_prev->velocity_x_, frame_prev->velocity_y_, frame_prev->velocity_z_, frame->velocity_x_, frame->velocity_y_);


		advect(frame->dt_, velocity_x, frame->velocity_x_, frame_prev->velocity_x_, frame_prev->velocity_x_, frame_prev->velocity_y_, frame_prev->velocity_z_);
		advect(frame->dt_, velocity_y, frame->velocity_y_, frame_prev->velocity_y_, frame_prev->velocity_x_, frame_prev->velocity_y_, frame_prev->velocity_z_);
		advect(frame->dt_, velocity_z, frame->velocity_z_, frame_prev->velocity_z_, frame_prev->velocity_x_, frame_prev->velocity_y_, frame_prev->velocity_z_);

		project(frame->velocity_x_, frame->velocity_y_, frame->velocity_z_, frame_prev->velocity_x_, frame_prev->velocity_y_);

		// density step

		diffuse(frame->dt_, frame->diff_, density, frame_prev->density_, frame->density_);
		advect(frame->dt_, density, frame->density_, frame_prev->density_, frame->velocity_x_, frame->velocity_y_, frame->velocity_z_);

		iter++;

		frame_prev->copy(n_, frame);
		//frame -> add dye/velocities
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

void FullAnimation::project(double* u, double* v, double* w, double* p, double* div) const
{
	unsigned long int i, j, k;
	//double* p = new double[52*52*52];
	//double* div = new double[52 * 52 * 52];
	double h = 1.0 / n_;

	for (i = 1; i <= n_; i++)
	{
		for (j = 1; j <= n_; j++)
		{
			for (k = 1; k <= n_; k++)
			{
				div[IDX(i, j, k, n_)] =
					-0.5f * h * (
						u[IDX(i + 1, j, k, n_)] -
						u[IDX(i - 1, j, k, n_)] +
						v[IDX(i, j + 1, k, n_)] -
						v[IDX(i, j - 1, k, n_)] +
						w[IDX(i, j, k + 1, n_)] -
						w[IDX(i, j, k - 1, n_)]);
				p[IDX(i, j, k, n_)] = 0;
			}
		}
	}
	// p - gradient field?
	set_bnd(density, div);
	set_bnd(density, p);
	lin_solve(density, p, div, 1, 6);
	
	for (i = 1; i <= n_; i++)
	{
		for (j = 1; j <= n_; j++)
		{
			for (k = 1; k <= n_; k++)
			{
				u[IDX(i, j, k, n_)] -= 0.5f * (
					p[IDX(i + 1, j, k, n_)] -
					p[IDX(i - 1, j, k, n_)]
					) * n_;
				v[IDX(i, j, k, n_)] -= 0.5f * (
					p[IDX(i, j + 1, k, n_)] -
					p[IDX(i, j - 1, k, n_)]
					) * n_;
				w[IDX(i, j, k, n_)] -= 0.5f * (
					p[IDX(i, j, k + 1, n_)] -
					p[IDX(i, j, k - 1, n_)]
					) * n_;
			}
		}
	}
	set_bnd(velocity_x, u);
	set_bnd(velocity_y, v);
	set_bnd(velocity_z, w);

	//delete[] div;
	//delete[] p;
}

void FullAnimation::set_bnd(int type, double* values) const
{
	unsigned long int i, j, k;
	for (i = 1; i <= n_; i++)
	{
		for (j = 1; j <= n_; j++)
		{
			values[IDX(i, j, 0, n_)] = (type == velocity_z ? -values[IDX(i, j, 1, n_)] : values[IDX(i, j, 1, n_)]);
			values[IDX(i, j, n_ + 1, n_)] = (type == velocity_z ? -values[IDX(i, j, n_, n_)] : values[IDX(i, j, n_, n_)]);
		}
	}
	for (i = 1; i <= n_; i++)
	{
		for (k = 1; k <= n_; k++)
		{
			values[IDX(i, 0, k, n_)] = (type == velocity_y ? -values[IDX(i, 1, k, n_)] : values[IDX(i, 1, k, n_)]);
			values[IDX(i, n_ + 1, k, n_)] = (type == velocity_y ? -values[IDX(i, n_, k, n_)] : values[IDX(i, n_, k, n_)]);
		}
	}
	for (j = 1; j <= n_; j++)
	{
		for (k = 1; k <= n_; k++)
		{
			values[IDX(0, j, k, n_)] = (type == velocity_x ? -values[IDX(1, j, k, n_)] : values[IDX(1, j, k, n_)]);
			values[IDX(n_ + 1, j, k, n_)] = (type == velocity_x ? -values[IDX(n_, j, k, n_)] : values[IDX(n_, j, k, n_)]);
		}
	}
	
	values[IDX(0, 0, 0, n_)] = 0.33 * (
		values[IDX(1, 0, 0, n_)] +
		values[IDX(0, 1, 0, n_)] +
		values[IDX(0, 0, 1, n_)]);
	values[IDX(0, n_ + 1, 0, n_)] = 0.33 * (
		values[IDX(1, n_ + 1, 0, n_)] +
		values[IDX(0, n_, 0, n_)] +
		values[IDX(0, n_ + 1, 1, n_)]);
	values[IDX(0, 0, n_ + 1, n_)] = 0.33 * (
		values[IDX(1, 0, n_ + 1, n_)] +
		values[IDX(0, 1, n_ + 1, n_)] +
		values[IDX(0, 0, n_, n_)]);
	values[IDX(0, n_ + 1, n_ + 1, n_)] = 0.33 * (
		values[IDX(1, n_ + 1, n_ + 1, n_)] +
		values[IDX(0, n_, n_ + 1, n_)] +
		values[IDX(0, n_ + 1, n_, n_)]);
	values[IDX(n_ + 1, 0, 0, n_)] = 0.33 * (
		values[IDX(n_, 0, 0, n_)] +
		values[IDX(n_ + 1, 1, 0, n_)] +
		values[IDX(n_ + 1, 0, 1, n_)]);
	values[IDX(n_ + 1, n_ + 1, 0, n_)] = 0.33 * (
		values[IDX(n_, n_ + 1, 0, n_)] +
		values[IDX(n_ + 1, n_, 0, n_)] +
		values[IDX(n_ + 1, n_ + 1, 1, n_)]);
	values[IDX(n_ + 1, 0, n_ + 1, n_)] = 0.33 * (
		values[IDX(n_, 0, n_ + 1, n_)] +
		values[IDX(n_ + 1, 1, n_ + 1, n_)] +
		values[IDX(n_ + 1, 0, n_, n_)]);
	values[IDX(n_ + 1, n_ + 1, n_ + 1, n_)] = 0.33 * (
		values[IDX(n_, n_ + 1, n_ + 1, n_)] +
		values[IDX(n_ + 1, n_, n_ + 1, n_)] +
		values[IDX(n_ + 1, n_ + 1, n_, n_)]);
}

void FullAnimation::diffuse(double dt, double diff, int type, double* values, double* values_0) const
{
	double a = dt * diff * n_ * n_ * n_;
	lin_solve(type, values, values_0, a, 1 + 6 * a);
}

void FullAnimation::advect(double dt, int type, double* d, double* d0, double* u, double* v, double* w)
{
	unsigned long int i, j, k;
	double x, y, z, s0, t0, r0, s1, t1, r1, dt0, i0, j0, k0, i1, j1, k1;

	dt0 = dt * n_;

	for (i = 1; i <= n_; i++)
	{
		for (j = 1; j <= n_; j++)
		{
			for (k = 1; k <= n_; k++)
			{
				x = static_cast<double>(i) - dt0 * u[IDX(i, j, k, n_)];
				y = static_cast<double>(j) - dt0 * v[IDX(i, j, k, n_)];
				z = static_cast<double>(k) - dt0 * w[IDX(i, j, k, n_)];
				if (x < 0.5) x = 0.5;
				if (x > n_ + 0.5) x = static_cast<double>(n_) + 0.5;
				i0 = floor(x);
				i1 = i0 + 1;
				if (y < 0.5) y = 0.5;
				if (y > n_ + 0.5) y = static_cast<double>(n_) + 0.5;
				j0 = floor(y);
				j1 = j0 + 1;
				if (z < 0.5) z = 0.5;
				if (z > n_ + 0.5) z = static_cast<double>(n_) + 0.5;
				k0 = floor(z);
				k1 = k0 + 1;
				
				s1 = x - i0;
				s0 = 1 - s1;
				t1 = y - j0;
				t0 = 1 - t1;
				r1 = z - k0;
				r0 = 1 - r1;

				int i0i = static_cast<int>(i0);
				int i1i = static_cast<int>(i1);
				int j0i = static_cast<int>(j0);
				int j1i = static_cast<int>(j1);
				int k0i = static_cast<int>(k0);
				int k1i = static_cast<int>(k1);

				const double res =
					s0 * t0 * r0 * d0[IDX(i0i, j0i, k0i, n_)] +
					s0 * t0 * r1 * d0[IDX(i0i, j0i, k1i, n_)] +
					s0 * t1 * r0 * d0[IDX(i0i, j1i, k0i, n_)] +
					s0 * t1 * r1 * d0[IDX(i0i, j1i, k1i, n_)] +
					s1 * t0 * r0 * d0[IDX(i1i, j0i, k0i, n_)] +
					s1 * t0 * r1 * d0[IDX(i1i, j0i, k1i, n_)] +
					s1 * t1 * r0 * d0[IDX(i1i, j1i, k0i, n_)] +
					s1 * t1 * r1 * d0[IDX(i1i, j1i, k1i, n_)];
									
				d[IDX(i, j, k, n_)] = res;

				if (type == density)
				{
					dens_volumes[iter*frame_size_+ IDX(i, j, k, n_)] = res; // TODO check order
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

	set_bnd(type, d);
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
					float dens_normal = static_cast<float>((dens_volumes[l*frame_size_ + IDX(i, j, k, n_)] - minTotal_) / (maxTotal_ - minTotal_));
					blender_file.write(reinterpret_cast<const char*>(&dens_normal), sizeof(dens_normal));
				}
			}
		}
	}

	
	blender_file.close();
}

FullAnimation::FullAnimation(): frame(nullptr), frame_prev(nullptr), dens_volumes(nullptr), time_(DEFAULT_time),
                                n_(DEFAULT_cube_side_size),
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
	delete frame;
	delete frame_prev;
	delete[] dens_volumes;
}
