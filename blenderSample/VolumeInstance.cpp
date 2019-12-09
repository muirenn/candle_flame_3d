#include "stdafx.h"
#include "VolumeInstance.h"


unsigned int VolumeInstance::idx(int i, int j, int k) const
{
	return (i + j * (n_ +2) + k * (n_ + 2) * (n_ + 2));
}

void VolumeInstance::write_to_file(float* arr)
{
	//float min = arr[idx(1,1,1)], max = arr[idx(1,1,1)];

	ofstream blender_file;
	blender_file.open(filename_, ios::out | ios::binary | ios::app);

	 find_minmax(arr);
	
	for (int i = 1; i <= n_; i++)
	{
		for (int j = 1; j <= n_; j++)
		{
			for (int k = 1; k <= n_; k++)
			{
				// C++ writes files in native format. On a standard x86-derived PC native is little endian
				  float currentV = static_cast<float>(arr[idx(i, j, k)] - min_) / (max_ - min_);
					//			float currentV = static_cast<float>(arr[idx(i, j, k)]);

				// zi = (xi - min(x))/(max(x)-min(x));

				// cout << currentV << endl;
				blender_file.write(reinterpret_cast<const char*>(&currentV), sizeof(currentV));
			}
			// cout << endl;
		}
		// cout << endl;
	}

	blender_file.close();
}

void VolumeInstance::add_source(float* x, float* s)
{
	for (unsigned int i = 0; i < size_; i++)
	{
		x[i] += dt_ * s[i];
	}
}

//TODO review the below
void VolumeInstance::project(float* u, float* v, float* w, float* p, float* div)
{
	int i, j, k, l;
	float h = 1.0f / n_;

//	float* p;
//	float* div;

//	p = new float[size_];
//	div = new float[size_];

	for (i = 1; i <= n_; i++)
	{
		for (j = 1; j <= n_; j++)
		{
			for (k = 1; k <= n_; k++)
			{
				div[idx(i, j, k)] =
					-0.5f * h * (
						u[idx(i + 1, j, k)] -
						u[idx(i - 1, j, k)] +
						v[idx(i, j + 1, k)] -
						v[idx(i, j - 1, k)] +
						w[idx(i, j, k + 1)] -
						w[idx(i, j, k - 1)]);
				p[idx(i, j, k)] = 0;
			}
		}
	}

	set_bnd(property::density, div);
	set_bnd(property::density, p);

	for (l = 0; l < iter; l++)
	{
		for (i = 1; i <= n_; i++)
		{
			for (j = 1; j <= n_; j++)
			{
				for (k = 1; k <= n_; k++)
				{
					p[idx(i, j, k)] = (
						div[idx(i, j, k)] +
						p[idx(i - 1, j, k)] +
						p[idx(i + 1, j, k)] +
						p[idx(i, j - 1, k)] +
						p[idx(i, j + 1, k)] +
						p[idx(i, j, k - 1)] +
						p[idx(i, j, k + 1)]
					) / 6;
				}
			}
		}
		set_bnd(property::density, p);
	}
	for (i = 1; i <= n_; i++)
	{
		for (j = 1; j <= n_; j++)
		{
			for (k = 1; k <= n_; k++)
			{
				u[idx(i, j, k)] -= 0.5f * (
					p[idx(i + 1, j, k)] -
					p[idx(i - 1, j, k)]
				) / h;
				v[idx(i, j, k)] -= 0.5f * (
					p[idx(i, j + 1, k)] -
					p[idx(i, j - 1, k)]
				) / h;
				w[idx(i, j, k)] -= 0.5f * (
					p[idx(i, j, k + 1)] -
					p[idx(i, j, k - 1)]
				) / h;
			}
		}
	}
	set_bnd(property::velocity_x, u);
	set_bnd(property::velocity_y, v);
	set_bnd(property::velocity_z, w);

	//delete[] p;
	//delete[] div;
}

void VolumeInstance::add_prev()
{
	int i, j, k;
	int max_dens = (n_ + 1) * (n_ + 1) * (n_ + 1);
	int half_n = (int)((n_ + 2) / 2);

	for (i = 0; i <= n_ + 1; i++)
	{
		for (j = 0; j <= n_ + 1; j++)
		{
			for (k = 0; k <= n_ + 1; k++)
			{
				float density_value = (i <= half_n ? i : n_ + 2 - i) * (j <= half_n ? j : n_+2 - j) * (k <= half_n ? k : n_+2 - k);
				density_prev_[idx(i, j, k)] += static_cast<float>(density_value);
				/*velocity_x_prev_[idx(i, j, k)] += avg_velocity;
				velocity_y_prev_[idx(i, j, k)] += high_velocity;
				velocity_z_prev_[idx(i, j, k)] += avg_velocity;*/
				
			}
		}
	}
}

void VolumeInstance::find_minmax(float* arr)
{
	// min_ = arr[0], max_ = arr[0];
	 min_ = 0, max_ = 4100;
	

	
	for (int i = 1; i <= n_; i++)
	{
		for (int j = 1; j <= n_; j++)
		{
			for (int k = 1; k <= n_; k++)
			{
				if (arr[idx(i, j, k)] > max_)
				{
					max_ = arr[idx(i, j, k)];
				}
				if (arr[idx(i, j, k)] < min_)
				{
					min_ = arr[idx(i, j, k)];
				}
			}
			// cout << endl;
		}
		// cout << endl;
	}
	
}

void VolumeInstance::set_bnd(int b, float* x) //TODO review
{
	/*
		int i;
		for (i = 1; i <= n_; i++) {
			x[IX(0, i)] = b == 1 ? –x[IX(1, i)] : x[IX(1, i)];
			x[IX(n_ + 1, i)] = b == 1 ? –x[IX(n_, i)] : x[IX(n_, i)];
			x[IX(i, 0)] = b == 2 ? –x[IX(i, 1)] : x[IX(i, 1)];
			x[IX(i, n_ + 1)] = b == 2 ? –x[IX(i, n_)] : x[IX(i, n_)];
	
			x[IX(0, 0, i)] = b == 1 ? –x[IX(0, 1, i)] : x[IX(0, 1, i)];
			x[IX(0, n_ + 1, i)] = b == 1 ? –x[IX(0, n_, i)] : x[IX(0, n_, i)];
			x[IX(n_ + 1, 0, i)] = b == 1 ? –x[IX(n_, 0, i)] : x[IX(n_, 0, i)];
			x[IX(n_ + 1, n_ + 1, i)] = b == 1 ? –x[IX(n_, i)] : x[IX(n_, i)];
	
			
			x[IX(i, 0)] = b == 2 ? –x[IX(i, 1)] : x[IX(i, 1)];
			x[IX(i, n_ + 1)] = b == 2 ? –x[IX(i, n_)] : x[IX(i, n_)];
		}
		*/
	
	for (int j = 1; j <= n_; j++)
	{
		for (int i = 1; i <= n_; i++)
		{
			x[idx(i, j, 0)] = b == property::velocity_z ? -x[idx(i, j, 1)] : x[idx(i, j, 1)];
			x[idx(i, j, n_ + 1)] = b == property::velocity_z ? -x[idx(i, j, n_)] : x[idx(i, j, n_)];
		}
	}
	for (int k = 1; k <= n_; k++)
	{
		for (int i = 1; i <= n_; i++)
		{
			x[idx(i, 0, k)] = b == property::velocity_y ? -x[idx(i, 1, k)] : x[idx(i, 1, k)];
			x[idx(i, n_ + 1, k)] = b == property::velocity_y ? -x[idx(i, n_, k)] : x[idx(i, n_, k)];
		}
	}
	for (int k = 1; k <= n_; k++)
	{
		for (int j = 1; j <= n_; j++)
		{
			x[idx(0, j, k)] = b == property::velocity_x ? -x[idx(1, j, k)] : x[idx(1, j, k)];
			x[idx(n_ + 1, j, k)] = b == property::velocity_x ? -x[idx(n_, j, k)] : x[idx(n_, j, k)];
		}
	}
	/*
	x[idx(0, 0, 0)] = (x[idx(1, 0, 0)] + x[idx(0, 1, 0)] + x[idx(0, 0, 1)]) / 3;
	x[idx(0, 0, n_ + 1)] = (x[idx(1, 0, n_ + 1)] + x[idx(0, 1, n_ + 1)] + x[idx(0, 0, n_)]) / 3;
	x[idx(0, n_ + 1, 0)] = (x[idx(1, n_ + 1, 0)] + x[idx(0, n_, 0)] + x[idx(0, n_ + 1, 1)]) / 3;
	x[idx(0, n_ + 1, n_ + 1)] = (x[idx(1, n_ + 1, n_ + 1)] + x[idx(0, n_, n_ + 1)] + x[idx(0, n_ + 1, n_)]) / 3;

	x[idx(n_ + 1, 0, 0)] = (x[idx(n_, 0, 0)] + x[idx(n_ + 1, 1, 0)] + x[idx(n_ + 1, 0, 1)]) / 3;
	x[idx(n_ + 1, 0, n_ + 1)] = (x[idx(n_, 0, n_ + 1)] + x[idx(n_ + 1, 1, n_ + 1)] + x[idx(n_ + 1, 0, n_)]) / 3;
	x[idx(n_ + 1, n_ + 1, 0)] = (x[idx(n_, n_ + 1, 0)] + x[idx(n_ + 1, n_, 0)] + x[idx(n_ + 1, n_ + 1, 1)]) / 3;
	x[idx(n_ + 1, n_ + 1, n_ + 1)] = (x[idx(n_, n_ + 1, n_ + 1)] + x[idx(n_ + 1, n_, n_ + 1)] + x[idx(
		n_ + 1, n_ + 1, n_)]) / 3;
	*/

	x[idx(0, 0, 0)] = 0.33f * (x[idx(1, 0, 0)]
		+ x[idx(0, 1, 0)]
		+ x[idx(0, 0, 1)]);
	x[idx(0, n_ + 1, 0)] = 0.33f * (x[idx(1, n_ + 1, 0)]
		+ x[idx(0, n_, 0)]
		+ x[idx(0, n_ + 1, 1)]);
	x[idx(0, 0, n_ + 1)] = 0.33f * (x[idx(1, 0, n_ + 1)]
		+ x[idx(0, 1, n_ + 1)]
		+ x[idx(0, 0, n_ + 2)]);
	x[idx(0, n_ + 1, n_ + 1)] = 0.33f * (x[idx(1, n_ + 1, n_ + 1)]
		+ x[idx(0, n_ + 2, n_ + 1)]
		+ x[idx(0, n_ + 1, n_ )]);
	x[idx(n_ + 1, 0, 0)] = 0.33f * (x[idx(n_, 0, 0)]
		+ x[idx(n_ + 1, 1, 0)]
		+ x[idx(n_ + 1, 0, 1)]);
	x[idx(n_ + 1, n_ + 1, 0)] = 0.33f * (x[idx(n_, n_ + 1, 0)]
		+ x[idx(n_ + 1, n_, 0)]
		+ x[idx(n_ + 1, n_ + 1, 1)]);
	x[idx(n_ + 1, 0, n_ + 1)] = 0.33f * (x[idx(n_, 0, n_ + 1)]
		+ x[idx(n_ + 1, 1, n_ + 1)]
		+ x[idx(n_ + 1, 0, n_)]);
	x[idx(n_ + 1, n_ + 1, n_ + 1)] = 0.33f * (x[idx(n_, n_ + 1, n_ + 1)]
		+ x[idx(n_ + 1, n_, n_ + 1)]
		+ x[idx(n_ + 1, n_ + 1, n_)]);
}

void VolumeInstance::diffuse(int b, float* x, float* x0)
{
	int i, j, k, l;
	float a = dt_ * diff_ * n_ * n_ * n_;
	for (l = 0; l < iter; l++)
	{
		for (i = 1; i <= n_; i++)
		{
			for (j = 1; j <= n_; j++)
			{
				for (k = 1; k <= n_; k++)
				{
					x[idx(i, j, k)] = (
						x0[idx(i - 1, j, k)] +
						x0[idx(i + 1, j, k)] +
						x0[idx(i, j - 1, k)] +
						x0[idx(i, j + 1, k)] +
						x0[idx(i, j, k - 1)] +
						x0[idx(i, j, k + 1)]
					) / (1 + 6 * a);
				}
			}
		}
		set_bnd(b, x);
	}
}

void VolumeInstance::advect(int b, float* d, float* d0, float* u, float* v, float* w)
{
	int i, j, k, i0, j0, k0, i1, j1, k1;
	float x, y, z, s0, t0, r0, s1, t1, r1, dt0;

	dt0 = dt_ * n_;

	for (i = 1; i <= n_; i++)
	{
		for (j = 1; j <= n_; j++)
		{
			for (k = 1; k <= n_; k++)
			{
				x = i - dt0 * u[idx(i, j, k)];
				y = j - dt0 * v[idx(i, j, k)];
				z = k - dt0 * w[idx(i, j, k)];
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

				d[idx(i, j, k)] =
					s0 * t0 * r0 * d0[idx(i0, j0, k0)] +
					s0 * t0 * r1 * d0[idx(i0, j0, k1)] +
					s0 * t1 * r0 * d0[idx(i0, j1, k0)] +
					s0 * t1 * r1 * d0[idx(i0, j1, k1)] +
					s1 * t0 * r0 * d0[idx(i1, j0, k0)] +
					s1 * t0 * r1 * d0[idx(i1, j0, k1)] +
					s1 * t1 * r0 * d0[idx(i1, j1, k0)] +
					s1 * t1 * r1 * d0[idx(i1, j1, k1)];
			}
		}
	}

	set_bnd(b, d);
}

void VolumeInstance::vel_step(float* velocity_x, float* velocity_y, float* velocity_z, float* velocity_x_prev,
                              float* velocity_y_prev, float* velocity_z_prev)
{
	add_source(velocity_x, velocity_x_prev);
	add_source(velocity_y, velocity_y_prev);
	add_source(velocity_z, velocity_z_prev);

	diffuse(1, velocity_x_prev, velocity_x);
	diffuse(2, velocity_y_prev, velocity_y);
	diffuse(3, velocity_z_prev, velocity_z);

	project(velocity_x_prev, velocity_y_prev, velocity_z_prev, velocity_x, velocity_y);

	advect(1, velocity_x, velocity_x_prev, velocity_x_prev, velocity_y_prev, velocity_z_prev);
	advect(2, velocity_y, velocity_y_prev, velocity_x_prev, velocity_y_prev, velocity_z_prev);
	advect(3, velocity_z, velocity_z_prev, velocity_x_prev, velocity_y_prev, velocity_z_prev);

	project(velocity_x, velocity_y, velocity_z, velocity_x_prev, velocity_y_prev);
}

void VolumeInstance::dens_step(float* density, float* density_prev, float* velocity_x, float* velocity_y,
                               float* velocity_z)
{
	add_source(density, density_prev);
	diffuse(0, density_prev, density);
	advect(0, density, density_prev, velocity_x, velocity_y, velocity_z);
}

VolumeInstance::VolumeInstance(): n_(cube_side_size), size_(0), visc_(viscosity), diff_(diffusivity), dt_(time_step),
                                  density_(nullptr),
                                  density_prev_(nullptr),
                                  velocity_x_(nullptr),
                                  velocity_y_(nullptr),
                                  velocity_z_(nullptr),
                                  velocity_x_prev_(nullptr),
                                  velocity_y_prev_(nullptr),
                                  velocity_z_prev_(nullptr), filename_(bvox_filename)
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

void VolumeInstance::init()
{
	size_ = (n_ + 2) * (n_ + 2) * (n_ + 2);

	density_ = new float[size_];
	density_prev_ = new float[size_];
	velocity_x_ = new float[size_]; //u in the script
	velocity_y_ = new float[size_]; //v
	velocity_z_ = new float[size_];
	velocity_x_prev_ = new float[size_];
	velocity_y_prev_ = new float[size_];
	velocity_z_prev_ = new float[size_];

	int i, j, k;
	int half_n = (int)((n_ + 2) / 2);
	int max_dens = (n_ + 1) * (n_ + 1) * (n_ + 1);

	char cc;
	for (i = 0; i <= n_+ 1; i++)
	{
		for (j = 0; j <= n_ + 1; j++)
		{
			for (k = 0; k <= n_ + 1; k++)
			{
				float density_value = (i <= half_n ? i : n_ + 2 - i) * (j <= half_n ? j : n_ + 2 - j) * (k <= half_n ? k : n_ + 2 - k); // center coloring trick
				// cout << "( " << i << " " << j << " " << k << " ) = " << density_value << endl;
			//	 cin >> cc;
				density_[idx(i, j, k)] = static_cast<float>(density_value);
				density_prev_[idx(i, j, k)] = static_cast<float>(density_value);
				velocity_x_[idx(i, j, k)] = avg_velocity;
				velocity_y_[idx(i, j, k)] = high_velocity;
				velocity_z_[idx(i, j, k)] = avg_velocity;
				velocity_x_prev_[idx(i, j, k)] = avg_velocity;
				velocity_y_prev_[idx(i, j, k)] = high_velocity;
				velocity_z_prev_[idx(i, j, k)] = avg_velocity;
			}
		}
	}
	/*
	const int strt = n_ % 9 * 4, fnsh = n_ % 9 * 6, avg = n_ % 9 * 5;
	 for (i = strt; i <= fnsh; i++)
	{
		for (j = strt; j <= fnsh; j++)
		{
			for (k = strt; k <= fnsh; k++)
			{
				density_prev_[idx(i, j, k)] = avg_density;
				density_[idx(i, j, k)] = avg_density;
			}
		}
	} 
	density_prev_[idx(avg, avg, avg)] = high_density;
	density_[idx(avg, avg, avg)] = high_density;
	*/

	ofstream blender_file;
	blender_file.open(filename_, ios::out | ios::binary | ios::trunc);

	blender_file.write(reinterpret_cast<const char*>(&n_), sizeof(n_));
	blender_file.write(reinterpret_cast<const char*>(&n_), sizeof(n_));
	blender_file.write(reinterpret_cast<const char*>(&n_), sizeof(n_));
	blender_file.close();
}

void VolumeInstance::step()
{
	// get_from_UI ( dens_prev, u_prev, v_prev );
	 vel_step(velocity_x_, velocity_y_, velocity_z_, velocity_x_prev_, velocity_y_prev_, velocity_z_prev_);
	 dens_step(density_, density_prev_, velocity_x_, velocity_y_, velocity_z_);
	write_to_file(density_);
	add_prev();
	// SWAP(velocity_x_prev_, velocity_x_);
	// SWAP(velocity_y_prev_, velocity_y_);
	// SWAP(velocity_z_prev_, velocity_z_);
	// SWAP(density_prev_, density_);
}
