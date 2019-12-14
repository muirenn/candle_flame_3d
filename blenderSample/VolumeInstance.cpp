#include "stdafx.h"
#include "VolumeInstance.h"

/*
unsigned int VolumeInstance::idx(int i, int j, int k) const
{
	return (i + j * (n_ +2) + k * (n_ + 2) * (n_ + 2));
}
*/
/*
void SWAP(double* x0, double* x)
{
//	double* tmp = x0;
//	x0 = x;
//	x = tmp;
}
*/

void VolumeInstance::write_to_file()
{
	//double min = arr[idx(1,1,1, n_)], max = arr[idx(1,1,1, n_)];

	ofstream blender_file;
	blender_file.open(filename_, ios::out | ios::binary | ios::app);

	// find_minmax(arr);
	min_ = 0; max_ = 5000;

	
	for (int i = 1; i <= n_; i++)
	{
		for (int j = 1; j <= n_; j++)
		{
			for (int k = 1; k <= n_; k++)
			{
				// C++ writes files in native format. On a standard x86-derived PC native is little endian
				  float currentV = static_cast<float>(density_[idx(i, j, k, n_)] - min_) / (max_ - min_);
				//				float currentV = static_cast<float>(arr[idx(i, j, k, n_)]);

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

void VolumeInstance::add_source(double* x, double* s)
{
	for (unsigned int i = 0; i < size_; i++)
	{
		x[i] += dt_ * s[i];
	}
}

//TODO review the below
void VolumeInstance::project(double* u, double* v, double* w, double* p, double* div)
{
	int i, j, k, l;
	double h = 1.0f / n_;

//	double* p;
//	double* div;

//	p = new double[size_];
//	div = new double[size_];

	for (i = 1; i <= n_; i++)
	{
		for (j = 1; j <= n_; j++)
		{
			for (k = 1; k <= n_; k++)
			{
				div[idx(i, j, k, n_)] =
					-0.5f * h * (
						u[idx(i + 1, j, k, n_)] -
						u[idx(i - 1, j, k, n_)] +
						v[idx(i, j + 1, k, n_)] -
						v[idx(i, j - 1, k, n_)] +
						w[idx(i, j, k + 1, n_)] -
						w[idx(i, j, k - 1, n_)]);
				p[idx(i, j, k, n_)] = 0;
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
					p[idx(i, j, k, n_)] = (
						div[idx(i, j, k, n_)] +
						p[idx(i - 1, j, k, n_)] +
						p[idx(i + 1, j, k, n_)] +
						p[idx(i, j - 1, k, n_)] +
						p[idx(i, j + 1, k, n_)] +
						p[idx(i, j, k - 1, n_)] +
						p[idx(i, j, k + 1, n_)]
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
				u[idx(i, j, k, n_)] -= 0.5f * (
					p[idx(i + 1, j, k, n_)] -
					p[idx(i - 1, j, k, n_)]
				) / h;
				v[idx(i, j, k, n_)] -= 0.5f * (
					p[idx(i, j + 1, k, n_)] -
					p[idx(i, j - 1, k, n_)]
				) / h;
				w[idx(i, j, k, n_)] -= 0.5f * (
					p[idx(i, j, k + 1, n_)] -
					p[idx(i, j, k - 1, n_)]
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
				//double density_value = (i <= half_n ? i : n_ + 2 - i) * (j <= half_n ? j : n_+2 - j) * (k <= half_n ? k : n_+2 - k);
				density_prev_[idx(i, j, k, n_)] = density_[idx(i, j, k, n_)]; // static_cast<double>(density_value);
				velocity_x_prev_[idx(i, j, k, n_)] = velocity_x_[idx(i, j, k, n_)];// avg_velocity;
				velocity_y_prev_[idx(i, j, k, n_)] = velocity_y_[idx(i, j, k, n_)]; // high_velocity;
				velocity_z_prev_[idx(i, j, k, n_)] = velocity_z_[idx(i, j, k, n_)];// avg_velocity;
				
			}
		}
	}
}

void VolumeInstance::find_minmax(double* arr)
{
	//
	 min_ = 0, max_ = 3e+9;
	

	/*
	min_ = arr[0], max_ = arr[0];
	for (int i = 1; i <= n_; i++)
	{
		for (int j = 1; j <= n_; j++)
		{
			for (int k = 1; k <= n_; k++)
			{
				if (arr[idx(i, j, k, n_)] > max_)
				{
					max_ = arr[idx(i, j, k, n_)];
				}
				if (arr[idx(i, j, k, n_)] < min_)
				{
					min_ = arr[idx(i, j, k, n_)];
				}
			}
			// cout << endl;
		}
		// cout << endl;
	}
*/	
}

void VolumeInstance::set_bnd(int b, double* x) //TODO review
{
	
	for (int j = 1; j <= n_; j++)
	{
		for (int i = 1; i <= n_; i++)
		{
			x[idx(i, j, 0, n_)] = b == property::velocity_z ? -x[idx(i, j, 1, n_)] : x[idx(i, j, 1, n_)];
			x[idx(i, j, n_ + 1, n_)] = b == property::velocity_z ? -x[idx(i, j, n_, n_)] : x[idx(i, j, n_, n_)];
		}
	}
	for (int k = 1; k <= n_; k++)
	{
		for (int i = 1; i <= n_; i++)
		{
			x[idx(i, 0, k, n_)] = b == property::velocity_y ? -x[idx(i, 1, k, n_)] : x[idx(i, 1, k, n_)];
			x[idx(i, n_ + 1, k, n_)] = b == property::velocity_y ? -x[idx(i, n_, k, n_)] : x[idx(i, n_, k, n_)];
		}
	}
	for (int k = 1; k <= n_; k++)
	{
		for (int j = 1; j <= n_; j++)
		{
			x[idx(0, j, k, n_)] = b == property::velocity_x ? -x[idx(1, j, k, n_)] : x[idx(1, j, k, n_)];
			x[idx(n_ + 1, j, k, n_)] = b == property::velocity_x ? -x[idx(n_, j, k, n_)] : x[idx(n_, j, k, n_)];
		}
	}

	x[idx(0, 0, 0, n_)] = 0.33f * (x[idx(1, 0, 0, n_)]
		+ x[idx(0, 1, 0, n_)]
		+ x[idx(0, 0, 1, n_)]);
	x[idx(0, n_ + 1, 0, n_)] = 0.33f * (x[idx(1, n_ + 1, 0, n_)]
		+ x[idx(0, n_, 0, n_)]
		+ x[idx(0, n_ + 1, 1, n_)]);
	x[idx(0, 0, n_ + 1, n_)] = 0.33f * (x[idx(1, 0, n_ + 1, n_)]
		+ x[idx(0, 1, n_ + 1, n_)]
		+ x[idx(0, 0, n_ + 2, n_)]);
	x[idx(0, n_ + 1, n_ + 1, n_)] = 0.33f * (x[idx(1, n_ + 1, n_ + 1, n_)]
		+ x[idx(0, n_ + 2, n_ + 1, n_)]
		+ x[idx(0, n_ + 1, n_ , n_)]);
	x[idx(n_ + 1, 0, 0, n_)] = 0.33f * (x[idx(n_, 0, 0, n_)]
		+ x[idx(n_ + 1, 1, 0, n_)]
		+ x[idx(n_ + 1, 0, 1, n_)]);
	x[idx(n_ + 1, n_ + 1, 0, n_)] = 0.33f * (x[idx(n_, n_ + 1, 0, n_)]
		+ x[idx(n_ + 1, n_, 0, n_)]
		+ x[idx(n_ + 1, n_ + 1, 1, n_)]);
	x[idx(n_ + 1, 0, n_ + 1, n_)] = 0.33f * (x[idx(n_, 0, n_ + 1, n_)]
		+ x[idx(n_ + 1, 1, n_ + 1, n_)]
		+ x[idx(n_ + 1, 0, n_, n_)]);
	x[idx(n_ + 1, n_ + 1, n_ + 1, n_)] = 0.33f * (x[idx(n_, n_ + 1, n_ + 1, n_)]
		+ x[idx(n_ + 1, n_, n_ + 1, n_)]
		+ x[idx(n_ + 1, n_ + 1, n_, n_)]);
}

void VolumeInstance::diffuse(int b, double* x, double* x0)
{
	int i, j, k, l;
	double a = dt_ * diff_ * n_ * n_ * n_;
	for (l = 0; l < iter; l++)
	{
		for (i = 1; i <= n_; i++)
		{
			for (j = 1; j <= n_; j++)
			{
				for (k = 1; k <= n_; k++)
				{
					x[idx(i, j, k, n_)] = x0[idx(i, j, k, n_)] + a*
						(
						x[idx(i - 1, j, k, n_)] +
						x[idx(i + 1, j, k, n_)] +
						x[idx(i, j - 1, k, n_)] +
						x[idx(i, j + 1, k, n_)] +
						x[idx(i, j, k - 1, n_)] +
						x[idx(i, j, k + 1, n_)]
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
				x = i - dt0 * u[idx(i, j, k, n_)];
				y = j - dt0 * v[idx(i, j, k, n_)];
				z = k - dt0 * w[idx(i, j, k, n_)];
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

				d[idx(i, j, k, n_)] =
					s0 * t0 * r0 * d0[idx(i0, j0, k0, n_)] +
					s0 * t0 * r1 * d0[idx(i0, j0, k1, n_)] +
					s0 * t1 * r0 * d0[idx(i0, j1, k0, n_)] +
					s0 * t1 * r1 * d0[idx(i0, j1, k1, n_)] +
					s1 * t0 * r0 * d0[idx(i1, j0, k0, n_)] +
					s1 * t0 * r1 * d0[idx(i1, j0, k1, n_)] +
					s1 * t1 * r0 * d0[idx(i1, j1, k0, n_)] +
					s1 * t1 * r1 * d0[idx(i1, j1, k1, n_)];
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

	density_ = new double[size_];
	density_prev_ = new double[size_];
	velocity_x_ = new double[size_]; //u in the script
	velocity_y_ = new double[size_]; //v
	velocity_z_ = new double[size_];
	velocity_x_prev_ = new double[size_];
	velocity_y_prev_ = new double[size_];
	velocity_z_prev_ = new double[size_];

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
				double density_value = (i <= half_n ? i : n_ + 2 - i) * (j <= half_n ? j : n_ + 2 - j) * (k <= half_n ? k : n_ + 2 - k); // center coloring trick
				// cout << "( " << i << " " << j << " " << k << " ) = " << density_value << endl;
			//	 cin >> cc;
				density_[idx(i, j, k, n_)] = static_cast<double>(density_value);
				density_prev_[idx(i, j, k, n_)] = static_cast<double>(density_value);
				velocity_x_[idx(i, j, k, n_)] = avg_velocity;
				velocity_y_[idx(i, j, k, n_)] = high_velocity;
				velocity_z_[idx(i, j, k, n_)] = avg_velocity;
				velocity_x_prev_[idx(i, j, k, n_)] = avg_velocity;
				velocity_y_prev_[idx(i, j, k, n_)] = high_velocity;
				velocity_z_prev_[idx(i, j, k, n_)] = avg_velocity;
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
				density_prev_[idx(i, j, k, n_)] = avg_density;
				density_[idx(i, j, k, n_)] = avg_density;
			}
		}
	} 
	density_prev_[idx(avg, avg, avg, n_)] = high_density;
	density_[idx(avg, avg, avg, n_)] = high_density;
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
	 vel_step();
	 dens_step();
	write_to_file();
	add_prev();
	// SWAP(velocity_x_prev_, velocity_x_);
	// SWAP(velocity_y_prev_, velocity_y_);
	// SWAP(velocity_z_prev_, velocity_z_);
	// SWAP(density_prev_, density_);
}
