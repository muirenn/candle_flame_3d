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
	const double h = 1.0f / n_;

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

	set_bnd(density, div);
	set_bnd(density, p);

	for (l = 0; l < DEFAULT_iter; l++)
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
	}
	for (i = 1; i <= n_; i++)
	{
		for (j = 1; j <= n_; j++)
		{
			for (k = 1; k <= n_; k++)
			{
				u[IDX(i, j, k, n_)] -= 0.5f * (
					p[IDX(i + 1, j, k, n_)] -
					p[IDX(i - 1, j, k, n_)]
				) / h;
				v[IDX(i, j, k, n_)] -= 0.5f * (
					p[IDX(i, j + 1, k, n_)] -
					p[IDX(i, j - 1, k, n_)]
				) / h;
				w[IDX(i, j, k, n_)] -= 0.5f * (
					p[IDX(i, j, k + 1, n_)] -
					p[IDX(i, j, k - 1, n_)]
				) / h;
			}
		}
	}
	set_bnd(velocity_x, u);
	set_bnd(velocity_y, v);
	set_bnd(velocity_z, w);
}

void VolumeInstance::set_bnd(int b, double* x) const
{
	for (int j = 1; j <= n_; j++)
	{
		for (int i = 1; i <= n_; i++)
		{
			x[IDX(i, j, 0, n_)] = b == velocity_z ? -x[IDX(i, j, 1, n_)] : x[IDX(i, j, 1, n_)];
			x[IDX(i, j, n_ + 1, n_)] = b == velocity_z ? -x[IDX(i, j, n_, n_)] : x[IDX(i, j, n_, n_)];
		}
	}
	for (int k = 1; k <= n_; k++)
	{
		for (int i = 1; i <= n_; i++)
		{
			x[IDX(i, 0, k, n_)] = b == velocity_y ? -x[IDX(i, 1, k, n_)] : x[IDX(i, 1, k, n_)];
			x[IDX(i, n_ + 1, k, n_)] = b == velocity_y ? -x[IDX(i, n_, k, n_)] : x[IDX(i, n_, k, n_)];
		}
	}
	for (int k = 1; k <= n_; k++)
	{
		for (int j = 1; j <= n_; j++)
		{
			x[IDX(0, j, k, n_)] = b == velocity_x ? -x[IDX(1, j, k, n_)] : x[IDX(1, j, k, n_)];
			x[IDX(n_ + 1, j, k, n_)] = b == velocity_x ? -x[IDX(n_, j, k, n_)] : x[IDX(n_, j, k, n_)];
		}
	}

	x[IDX(0, 0, 0, n_)] = 0.33f * (x[IDX(1, 0, 0, n_)]
		+ x[IDX(0, 1, 0, n_)]
		+ x[IDX(0, 0, 1, n_)]);
	x[IDX(0, n_ + 1, 0, n_)] = 0.33f * (x[IDX(1, n_ + 1, 0, n_)]
		+ x[IDX(0, n_, 0, n_)]
		+ x[IDX(0, n_ + 1, 1, n_)]);
	x[IDX(0, 0, n_ + 1, n_)] = 0.33f * (x[IDX(1, 0, n_ + 1, n_)]
		+ x[IDX(0, 1, n_ + 1, n_)]
		+ x[IDX(0, 0, n_ + 2, n_)]);
	x[IDX(0, n_ + 1, n_ + 1, n_)] = 0.33f * (x[IDX(1, n_ + 1, n_ + 1, n_)]
		+ x[IDX(0, n_ + 2, n_ + 1, n_)]
		+ x[IDX(0, n_ + 1, n_, n_)]);
	x[IDX(n_ + 1, 0, 0, n_)] = 0.33f * (x[IDX(n_, 0, 0, n_)]
		+ x[IDX(n_ + 1, 1, 0, n_)]
		+ x[IDX(n_ + 1, 0, 1, n_)]);
	x[IDX(n_ + 1, n_ + 1, 0, n_)] = 0.33f * (x[IDX(n_, n_ + 1, 0, n_)]
		+ x[IDX(n_ + 1, n_, 0, n_)]
		+ x[IDX(n_ + 1, n_ + 1, 1, n_)]);
	x[IDX(n_ + 1, 0, n_ + 1, n_)] = 0.33f * (x[IDX(n_, 0, n_ + 1, n_)]
		+ x[IDX(n_ + 1, 1, n_ + 1, n_)]
		+ x[IDX(n_ + 1, 0, n_, n_)]);
	x[IDX(n_ + 1, n_ + 1, n_ + 1, n_)] = 0.33f * (x[IDX(n_, n_ + 1, n_ + 1, n_)]
		+ x[IDX(n_ + 1, n_, n_ + 1, n_)]
		+ x[IDX(n_ + 1, n_ + 1, n_, n_)]);
}

void VolumeInstance::diffuse(int b, double* x, double* x0)
{
	int i, j, k, l;
	double a = dt_ * diff_ * n_ * n_ * n_;
	for (l = 0; l < DEFAULT_iter; l++)
	{
		for (i = 1; i <= n_; i++)
		{
			for (j = 1; j <= n_; j++)
			{
				for (k = 1; k <= n_; k++)
				{
					x[IDX(i, j, k, n_)] = x0[IDX(i, j, k, n_)] + a *
					(
						x[IDX(i - 1, j, k, n_)] +
						x[IDX(i + 1, j, k, n_)] +
						x[IDX(i, j - 1, k, n_)] +
						x[IDX(i, j + 1, k, n_)] +
						x[IDX(i, j, k - 1, n_)] +
						x[IDX(i, j, k + 1, n_)]
					) / (1 + 6 * a);
				}
			}
		}
		set_bnd(b, x);
	}
}

void VolumeInstance::advect(int b, double* d, double* d0, double* u, double* v, double* w)
{
	int i, j, k, i0, j0, k0, i1, j1, k1;
	double x, y, z, s0, t0, r0, s1, t1, r1, dt0;

	dt0 = dt_ * n_;

	for (i = 1; i <= n_; i++)
	{
		for (j = 1; j <= n_; j++)
		{
			for (k = 1; k <= n_; k++)
			{
				x = i - dt0 * u[IDX(i, j, k, n_)];
				y = j - dt0 * v[IDX(i, j, k, n_)];
				z = k - dt0 * w[IDX(i, j, k, n_)];
				if (x < 0.5) x = 0.5;
				if (x > n_ + 0.5) x = n_ + 0.5f;
				i0 = static_cast<int>(x);
				i1 = i0 + 1;
				if (y < 0.5) y = 0.5;
				if (y > n_ + 0.5) y = n_ + 0.5f;
				j0 = static_cast<int>(y);
				j1 = j0 + 1;
				if (z < 0.5) z = 0.5;
				if (z > n_ + 0.5) z = n_ + 0.5f;
				k0 = static_cast<int>(z);
				k1 = k0 + 1;
				s1 = x - i0;
				s0 = 1 - s1;
				t1 = y - j0;
				t0 = 1 - t1;
				r1 = z - k0;
				r0 = 1 - r1;

				const float dens =
					s0 * t0 * r0 * d0[IDX(i0, j0, k0, n_)] +
					s0 * t0 * r1 * d0[IDX(i0, j0, k1, n_)] +
					s0 * t1 * r0 * d0[IDX(i0, j1, k0, n_)] +
					s0 * t1 * r1 * d0[IDX(i0, j1, k1, n_)] +
					s1 * t0 * r0 * d0[IDX(i1, j0, k0, n_)] +
					s1 * t0 * r1 * d0[IDX(i1, j0, k1, n_)] +
					s1 * t1 * r0 * d0[IDX(i1, j1, k0, n_)] +
					s1 * t1 * r1 * d0[IDX(i1, j1, k1, n_)];

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

void VolumeInstance::vel_step()
{
	//	add_source(velocity_x_, velocity_x_prev_);
	//	add_source(velocity_y_, velocity_y_prev_);
	//add_source(velocity_z_, velocity_z_prev_);

	diffuse(1, velocity_x_prev_, velocity_x_);
	diffuse(2, velocity_y_prev_, velocity_y_);
	diffuse(3, velocity_z_prev_, velocity_z_);

	project(velocity_x_prev_, velocity_y_prev_, velocity_z_prev_, velocity_x_, velocity_y_);


	advect(1, velocity_x_, velocity_x_prev_, velocity_x_prev_, velocity_y_prev_, velocity_z_prev_);
	advect(2, velocity_y_, velocity_y_prev_, velocity_x_prev_, velocity_y_prev_, velocity_z_prev_);
	advect(3, velocity_z_, velocity_z_prev_, velocity_x_prev_, velocity_y_prev_, velocity_z_prev_);

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
	delete[] velocity_x_;
	velocity_x_ = nullptr;
	delete[] velocity_y_;
	velocity_y_ = nullptr;
	delete[] velocity_z_;
	velocity_z_ = nullptr;
	delete[] velocity_x_prev_;
	velocity_x_prev_ = nullptr;
	delete[] velocity_y_prev_;
	velocity_y_prev_ = nullptr;
	delete[] velocity_z_prev_;
	velocity_z_prev_ = nullptr;
	cout << endl << "Destructor called" << endl;
}

void VolumeInstance::allocate_memory()
{
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
	const int half_n = static_cast<int>((n_ + 2) / 2);

	for (i = 0; i <= n_ + 1; i++)
	{
		for (j = 0; j <= n_ + 1; j++)
		{
			for (k = 0; k <= n_ + 1; k++)
			{
				const double density_value = (i <= half_n ? i : n_ + 2 - i) * (j <= half_n ? j : n_ + 2 - j) * (
					k <= half_n ? k : n_ + 2 - k); // center coloring trick
				//density_[IDX(i, j, k, n_)] = static_cast<double>(density_value);
				density_prev_[IDX(i, j, k, n_)] = static_cast<double>(density_value);
				//velocity_x_[IDX(i, j, k, n_)] = DEFAULT_avg_velocity;
				//velocity_y_[IDX(i, j, k, n_)] = DEFAULT_high_velocity;
				//velocity_z_[IDX(i, j, k, n_)] = DEFAULT_avg_velocity;
				velocity_x_prev_[IDX(i, j, k, n_)] = DEFAULT_avg_velocity;
				velocity_y_prev_[IDX(i, j, k, n_)] = DEFAULT_high_velocity;
				velocity_z_prev_[IDX(i, j, k, n_)] = DEFAULT_avg_velocity;
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

void VolumeInstance::set_prev_values(const VolumeInstance& prev_frame)
{
	density_prev_ = prev_frame.density_;
	velocity_x_prev_ = prev_frame.velocity_x_;
	velocity_y_prev_ = prev_frame.velocity_y_;
	velocity_z_prev_ = prev_frame.velocity_z_;

	n_ = prev_frame.n_;
	size_ = prev_frame.size_;
	visc_ = prev_frame.visc_;
	diff_ = prev_frame.diff_;
	dt_ = prev_frame.dt_;
}
