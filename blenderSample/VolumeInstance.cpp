#include "stdafx.h"
#include "VolumeInstance.h"

VolumeInstance::VolumeInstance(): visc_(DEFAULT_viscosity),
                                  diff_(DEFAULT_diffusivity), dt_(DEFAULT_time_step), min_(0), max_(0),
                                  density_(nullptr),
                                  velocity_x_(nullptr),
                                  velocity_y_(nullptr),
                                  velocity_z_(nullptr)
{
}

VolumeInstance::~VolumeInstance()
{
	delete[] density_;
	density_ = nullptr;
	delete[] velocity_x_;
	velocity_x_ = nullptr;
	delete[] velocity_y_;
	velocity_y_ = nullptr;
	delete[] velocity_z_;
	velocity_z_ = nullptr;
}

void VolumeInstance::allocate_memory(unsigned long long int size)
{
	density_ = new double[size];
	velocity_x_ = new double[size];
	velocity_y_ = new double[size];
	velocity_z_ = new double[size];

	visc_ = DEFAULT_viscosity;
	diff_ = DEFAULT_diffusivity;
	dt_ = DEFAULT_time_step;
}

void VolumeInstance::draw_sphere(unsigned long int n)
{
	unsigned long int i, j, k;
	const unsigned long int half_n = static_cast<int>(n / 2);
	max_ = 0; min_ = 0;

	for (i = 0; i <= n + 1; i++)
	{
		for (j = 0; j <= n + 1; j++)
		{
			for (k = 0; k <= n + 1; k++)
			{
				const double density_value = static_cast<double>(i <= half_n ? i : n - i) * static_cast<double>(j <= half_n ? j : n - j) * static_cast<double>(
					k <= half_n ? k : n - k); 
				density_[IDX(i, j, k, n)] = static_cast<double>(density_value);
				velocity_x_[IDX(i, j, k, n)] = DEFAULT_high_velocity;
				velocity_y_[IDX(i, j, k, n)] = 0;
				velocity_z_[IDX(i, j, k, n)] = 0;

				if(density_value > max_)
				{
					max_ = density_value;
				}
				if (density_value < min_)
				{
					min_ = density_value;
				}
				
			}
		}
	}
}

void VolumeInstance::draw_candle_v1(unsigned long int n) const
{
	unsigned long int i, j, k;
	const unsigned long int half_n = static_cast<int>(n / 2);
	
	for (i = 0; i <= n + 1; i++)
	{
		for (j = 0; j <= n + 1; j++)
		{
			for (k = 0; k <= n + 1; k++)
			{
			//	density_[IDX(i, j, k, n_)] = 0;
				// density_prev_[IDX(i, j, k, n_)] = 0;
				velocity_x_[IDX(i, j, k, n)] = 0;
				velocity_y_[IDX(i, j, k, n)] = DEFAULT_avg_velocity;
				velocity_z_[IDX(i, j, k, n)] = 0;



				if(i == 0 || j == 0 || k == 0 || i == n + 1 || j == n +1 || k == n + 1)
				{
					density_[IDX(i, j, k, n)] = 0;
				}
				else
				{
					const int x_displaced = k - half_n - 1;
					const int y_displaced = j - half_n - 1;
					int rad_coord = static_cast<int>(sqrt(x_displaced * x_displaced + y_displaced * y_displaced)) - 1;
					if (rad_coord >= half_n)
					{
						density_[IDX(i, j, k, n)] = 0;
					}
					else
					{
						const int h_coord = n - i;
						density_[IDX(i, j, k, n)] = static_cast<double>(candle_80x80[h_coord * half_n + rad_coord]);
					}
				}
			}
		}
	}
}

double VolumeInstance::get_min() const
{
	return min_;
}

double VolumeInstance::get_max() const
{
	return max_;
}

void VolumeInstance::fill_with_zeros(unsigned long int n)
{
	unsigned long int i, j, k;
	max_ = 0; min_ = 0;

	for (i = 0; i <= n + 1; i++)
	{
		for (j = 0; j <= n + 1; j++)
		{
			for (k = 0; k <= n + 1; k++)
			{
				density_[IDX(i, j, k, n)] = 0;
				velocity_x_[IDX(i, j, k, n)] = 0;
				velocity_y_[IDX(i, j, k, n)] = 0;
				velocity_z_[IDX(i, j, k, n)] = 0;
			}
		}
	}
}
