#include "stdafx.h"
#include "VolumeInstance.h"

void VolumeInstance::add_source(double* x, double* s) const
{
	for (unsigned int i = 0; i < size_; i++)
	{
		x[i] += dt_ * s[i];
	}
}

void VolumeInstance::project(double* u, double* v, double* w, double* p, double* div) const
{
	int i, j, k, l;
	// const double h = 1.0f / n_;

	for (i = 1; i <= n_; i++)
	{
		for (j = 1; j <= n_; j++)
		{
			for (k = 1; k <= n_; k++)
			{
				div[IDX(i, j, k, n_)] =
					-0.5f / static_cast<double>(n_) * (
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

	set_bnd(density, div);
	set_bnd(density, p);

	lin_solve(density, p, div, 1, 6);
	/*for (l = 0; l < DEFAULT_iter; l++)
	{
		for (i = 1; i <= n_; i++)
		{
			for (j = 1; j <= n_; j++)
			{
				for (k = 1; k <= n_; k++)
				{
					p[IDX(i, j, k, n_)] = (
						div[IDX(i, j, k, n_)] +
						p[IDX(i - 1, j, k, n_)] +
						p[IDX(i + 1, j, k, n_)] +
						p[IDX(i, j - 1, k, n_)] +
						p[IDX(i, j + 1, k, n_)] +
						p[IDX(i, j, k - 1, n_)] +
						p[IDX(i, j, k + 1, n_)]
					) / 6;
				}
			}
		}
		set_bnd(density, p);
	}*/
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
}

void VolumeInstance::set_bnd(const int type, double* values) const
{
	for (int i = 1; i <= n_; i++)
	{
		for (int j = 1; j <= n_; j++)
		{
			values[IDX(i, j, 0, n_)] = (type == velocity_z ? -values[IDX(i, j, 1, n_)] : values[IDX(i, j, 1, n_)]);
			values[IDX(i, j, n_ + 1, n_)] = (type == velocity_z ? -values[IDX(i, j, n_, n_)] : values[IDX(i, j, n_, n_)]);
		}
	}
	for (int i = 1; i <= n_; i++)
	{
		for (int k = 1; k <= n_; k++)
		{
			values[IDX(i, 0, k, n_)] = (type == velocity_y ? -values[IDX(i, 1, k, n_)] : values[IDX(i, 1, k, n_)]);
			values[IDX(i, n_ + 1, k, n_)] = (type == velocity_y ? -values[IDX(i, n_, k, n_)] : values[IDX(i, n_, k, n_)]);
		}
	}
	for (int j = 1; j <= n_; j++)
	{
		for (int k = 1; k <= n_; k++)
		{
			values[IDX(0, j, k, n_)] = (type == velocity_x ? -values[IDX(1, j, k, n_)] : values[IDX(1, j, k, n_)]);
			values[IDX(n_ + 1, j, k, n_)] = (type == velocity_x ? -values[IDX(n_, j, k, n_)] : values[IDX(n_, j, k, n_)]);
		}
	}

	values[IDX(0, 0, 0, n_)] = 0.33f * (values[IDX(1, 0, 0, n_)]
		+ values[IDX(0, 1, 0, n_)]
		+ values[IDX(0, 0, 1, n_)]);
	values[IDX(0, n_ + 1, 0, n_)] = 0.33f * (values[IDX(1, n_ + 1, 0, n_)]
		+ values[IDX(0, n_, 0, n_)]
		+ values[IDX(0, n_ + 1, 1, n_)]);
	values[IDX(0, 0, n_ + 1, n_)] = 0.33f * (values[IDX(1, 0, n_ + 1, n_)]
		+ values[IDX(0, 1, n_ + 1, n_)]
		+ values[IDX(0, 0, n_ + 2, n_)]);
	values[IDX(0, n_ + 1, n_ + 1, n_)] = 0.33f * (values[IDX(1, n_ + 1, n_ + 1, n_)]
		+ values[IDX(0, n_ + 2, n_ + 1, n_)]
		+ values[IDX(0, n_ + 1, n_, n_)]);
	values[IDX(n_ + 1, 0, 0, n_)] = 0.33f * (values[IDX(n_, 0, 0, n_)]
		+ values[IDX(n_ + 1, 1, 0, n_)]
		+ values[IDX(n_ + 1, 0, 1, n_)]);
	values[IDX(n_ + 1, n_ + 1, 0, n_)] = 0.33f * (values[IDX(n_, n_ + 1, 0, n_)]
		+ values[IDX(n_ + 1, n_, 0, n_)]
		+ values[IDX(n_ + 1, n_ + 1, 1, n_)]);
	values[IDX(n_ + 1, 0, n_ + 1, n_)] = 0.33f * (values[IDX(n_, 0, n_ + 1, n_)]
		+ values[IDX(n_ + 1, 1, n_ + 1, n_)]
		+ values[IDX(n_ + 1, 0, n_, n_)]);
	values[IDX(n_ + 1, n_ + 1, n_ + 1, n_)] = 0.33f * (values[IDX(n_, n_ + 1, n_ + 1, n_)]
		+ values[IDX(n_ + 1, n_, n_ + 1, n_)]
		+ values[IDX(n_ + 1, n_ + 1, n_, n_)]);
}

void VolumeInstance::diffuse(const int type, double* values, double* values_0) const
{
	double a = dt_ * diff_ * n_ * n_ * n_;

	lin_solve(type, values, values_0, a, 1 + 6 * a);
	
	/*for (l = 0; l < DEFAULT_iter; l++)
	{
		for (i = 1; i <= n_; i++)
		{
			for (j = 1; j <= n_; j++)
			{
				for (k = 1; k <= n_; k++)
				{
					values[IDX(i, j, k, n_)] = values_0[IDX(i, j, k, n_)] + a *
					(
						values[IDX(i - 1, j, k, n_)] +
						values[IDX(i + 1, j, k, n_)] +
						values[IDX(i, j - 1, k, n_)] +
						values[IDX(i, j + 1, k, n_)] +
						values[IDX(i, j, k - 1, n_)] +
						values[IDX(i, j, k + 1, n_)]
					) / (1 + 6 * a);
				}
			}
		}
		set_bnd(type, values);
	}*/
}

void VolumeInstance::advect(int b, double* d, double* d0, double* u, double* v, double* w)
{
	int i, j, k;
	double x, y, z, s0, t0, r0, s1, t1, r1, dt0, i0, j0, k0, i1, j1, k1;

	dt0 = dt_ * n_;

	for (k = 1; k <= n_; k++)
	{
		for (j = 1; j <= n_; j++)
		{
			for (i = 1; i <= n_; i++)
			{
				x = static_cast<double>(i) - dt0 * u[IDX(i, j, k, n_)];
				y = static_cast<double>(j) - dt0 * v[IDX(i, j, k, n_)];
				z = static_cast<double>(k) - dt0 * w[IDX(i, j, k, n_)];
				if (x < 0.5) x = 0.5;
				if (x > n_ + 0.5) x = static_cast<double>(n_) + 0.5;
				i0 = floorf(x);
				i1 = i0 + 1;
				if (y < 0.5) y = 0.5;
				if (y > n_ + 0.5) y = static_cast<double>(n_) + 0.5;
				j0 = floorf(y);
				j1 = j0 + 1;
				if (z < 0.5) z = 0.5;
				if (z > n_ + 0.5) z = static_cast<double>(n_) + 0.5;
				k0 = floorf(z);
				k1 = k0 + 1;
				s1 = x - i0;
				s0 = 1 - s1;
				t1 = y - j0;
				t0 = 1 - t1;
				r1 = z - k0;
				r0 = 1 - r1;

				int i0i = i0;
				int i1i = i1;
				int j0i = j0;
				int j1i = j1;
				int k0i = k0;
				int k1i = k1;

				const double dens =
					s0 * t0 * r0 * d0[IDX(i0i, j0i, k0i, n_)] +
					s0 * t0 * r1 * d0[IDX(i0i, j0i, k1i, n_)] +
					s0 * t1 * r0 * d0[IDX(i0i, j1i, k0i, n_)] +
					s0 * t1 * r1 * d0[IDX(i0i, j1i, k1i, n_)] +
					s1 * t0 * r0 * d0[IDX(i1i, j0i, k0i, n_)] +
					s1 * t0 * r1 * d0[IDX(i1i, j0i, k1i, n_)] +
					s1 * t1 * r0 * d0[IDX(i1i, j1i, k0i, n_)] +
					s1 * t1 * r1 * d0[IDX(i1i, j1i, k1i, n_)];

				d[IDX(i, j, k, n_)] = dens;


				if (b == density)
				{
					if (dens < min_)
					{
						min_ = dens;
					}
					if (dens > max_)
					{
						max_ = dens;
					}
				}
			}
		}
	}

	set_bnd(b, d);
}

void VolumeInstance::lin_solve(int type, double* values, double* values_0, double a, double c) const
{
	int i, j, k, l;
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
							)) / c;
				}
			}
		}
		set_bnd(type, values);
	}
}

void VolumeInstance::vel_step()
{
	//	add_source(velocity_x_, velocity_x_prev_);
	//	add_source(velocity_y_, velocity_y_prev_);
	//add_source(velocity_z_, velocity_z_prev_);

	diffuse(velocity_x, velocity_x_prev_, velocity_x_);
	diffuse(velocity_y, velocity_y_prev_, velocity_y_);
	diffuse(velocity_z, velocity_z_prev_, velocity_z_);

	project(velocity_x_prev_, velocity_y_prev_, velocity_z_prev_, velocity_x_, velocity_y_);


	advect(velocity_x, velocity_x_, velocity_x_prev_, velocity_x_prev_, velocity_y_prev_, velocity_z_prev_);
	advect(velocity_y, velocity_y_, velocity_y_prev_, velocity_x_prev_, velocity_y_prev_, velocity_z_prev_);
	advect(velocity_z, velocity_z_, velocity_z_prev_, velocity_x_prev_, velocity_y_prev_, velocity_z_prev_);

	project(velocity_x_, velocity_y_, velocity_z_, velocity_x_prev_, velocity_y_prev_);
}

void VolumeInstance::dens_step()
{
	//add_source(density_, density_prev_);
	diffuse(0, density_prev_, density_);
	advect(0, density_, density_prev_, velocity_x_, velocity_y_, velocity_z_);
}

VolumeInstance::VolumeInstance(): n_(DEFAULT_cube_side_size), size_(0), visc_(DEFAULT_viscosity),
                                  diff_(DEFAULT_diffusivity), dt_(DEFAULT_time_step), min_(0), max_(0),
                                  density_(nullptr),
                                  density_prev_(nullptr),
                                  velocity_x_(nullptr),
                                  velocity_y_(nullptr),
                                  velocity_z_(nullptr),
                                  velocity_x_prev_(nullptr),
                                  velocity_y_prev_(nullptr),
                                  velocity_z_prev_(nullptr)
{
}

VolumeInstance::~VolumeInstance()
{
	delete[] density_;
	density_ = nullptr;
		delete[] density_prev_;
		density_prev_ = nullptr;
		delete[] velocity_x_prev_;
		velocity_x_prev_ = nullptr;
		delete[] velocity_y_prev_;
		velocity_y_prev_ = nullptr;
		delete[] velocity_z_prev_;
		velocity_z_prev_ = nullptr;

	
	delete[] velocity_x_;
	velocity_x_ = nullptr;
	delete[] velocity_y_;
	velocity_y_ = nullptr;
	delete[] velocity_z_;
	velocity_z_ = nullptr;
}

void VolumeInstance::allocate_memory(int n)
{
	n_ = n;
	size_ = (n_ + 2) * (n_ + 2) * (n_ + 2);

	density_ = new double[size_];
	density_prev_ = new double[size_];
	velocity_x_ = new double[size_];
	velocity_y_ = new double[size_];
	velocity_z_ = new double[size_];
	velocity_x_prev_ = new double[size_];
	velocity_y_prev_ = new double[size_];
	velocity_z_prev_ = new double[size_];
}

void VolumeInstance::step()
{
	vel_step();
	dens_step();
}

void VolumeInstance::draw_sphere() const
{
	int i, j, k;
	const int half_n = static_cast<int>(n_ / 2);

	for (i = 0; i <= n_ + 1; i++)
	{
		for (j = 0; j <= n_ + 1; j++)
		{
			for (k = 0; k <= n_ + 1; k++)
			{
				const double density_value = static_cast<double>(i <= half_n ? i : n_ - i) * static_cast<double>(j <= half_n ? j : n_ - j) * static_cast<double>(
					k <= half_n ? k : n_ - k); 
				density_[IDX(i, j, k, n_)] = static_cast<double>(density_value);
				density_prev_[IDX(i, j, k, n_)] = static_cast<double>(density_value);
				velocity_x_[IDX(i, j, k, n_)] = DEFAULT_avg_velocity;
				velocity_y_[IDX(i, j, k, n_)] = DEFAULT_high_velocity;
				velocity_z_[IDX(i, j, k, n_)] = DEFAULT_avg_velocity;
				velocity_x_prev_[IDX(i, j, k, n_)] = DEFAULT_avg_velocity;
				velocity_y_prev_[IDX(i, j, k, n_)] = DEFAULT_high_velocity;
				velocity_z_prev_[IDX(i, j, k, n_)] = DEFAULT_avg_velocity;
			}
		}
	}
}

void VolumeInstance::draw_candle_v1() const
{
	int i, j, k;
	const int half_n = static_cast<int>(n_ / 2);
	
	for (i = 0; i <= n_ + 1; i++)
	{
		for (j = 0; j <= n_ + 1; j++)
		{
			for (k = 0; k <= n_ + 1; k++)
			{
			//	density_[IDX(i, j, k, n_)] = 0;
				// density_prev_[IDX(i, j, k, n_)] = 0;
				velocity_x_[IDX(i, j, k, n_)] = DEFAULT_avg_velocity;
				velocity_y_[IDX(i, j, k, n_)] = 0;
				velocity_z_[IDX(i, j, k, n_)] = 0;
				velocity_x_prev_[IDX(i, j, k, n_)] = DEFAULT_avg_velocity;
				velocity_y_prev_[IDX(i, j, k, n_)] = 0;
				velocity_z_prev_[IDX(i, j, k, n_)] = 0;


				if(i == 0 || j == 0 || k == 0 || i == n_ + 1 || j == n_ +1 || k == n_ + 1)
				{
					density_[IDX(i, j, k, n_)] = 0;
				}
				else
				{
					const int x_displaced = k - half_n - 1;
					const int y_displaced = j - half_n - 1;
					int rad_coord = static_cast<int>(sqrt(x_displaced * x_displaced + y_displaced * y_displaced)) - 1;
					if (rad_coord >= half_n)
					{
						density_[IDX(i, j, k, n_)] = 0;
					}
					else
					{
						const int h_coord = n_ - i;
						density_[IDX(i, j, k, n_)] = static_cast<double>(candle_80x80[h_coord * half_n + rad_coord]);
					}

					//cout << density_[IDX(i, j, k, n_)] << endl;
				}
				density_prev_[IDX(i, j, k, n_)] = density_[IDX(i, j, k, n_)];

			}
		}
	}
}

double VolumeInstance::get_density(const int i) const
{
	return density_[i];
}

double VolumeInstance::get_min() const
{
	return min_;
}

double VolumeInstance::get_max() const
{
	return max_;
}

void VolumeInstance::set_prev_values_nochange(const VolumeInstance& prev_frame)
{
	int i, j, k;
	for (i = 0; i <= n_ + 1; i++)
	{
		for (j = 0; j <= n_ + 1; j++)
		{
			for (k = 0; k <= n_ + 1; k++)
			{
				int idx = IDX(i, j, k, n_);
				density_prev_[idx] = prev_frame.density_[idx];
				velocity_x_prev_[idx] = prev_frame.velocity_x_[idx];
				velocity_y_prev_[idx] = prev_frame.velocity_y_[idx];
				velocity_z_prev_[idx] = prev_frame.velocity_z_[idx];
				density_[idx] = prev_frame.density_[idx];
				velocity_x_[idx] = prev_frame.velocity_x_[idx];
				velocity_y_[idx] = prev_frame.velocity_y_[idx];
				velocity_z_[idx] = prev_frame.velocity_z_[idx];
			}
		}
	}
	


	n_ = prev_frame.n_;
	size_ = prev_frame.size_;
	visc_ = prev_frame.visc_;
	diff_ = prev_frame.diff_;
	dt_ = prev_frame.dt_;
}
