#include "stdafx.h"
#include "VolumeInstance.h"

VolumeInstance::VolumeInstance(): visc_(AIR_VISC),
                                  diff_(AIR_DIFF), dt_(DEFAULT_time_step), min_(0), max_(0),
                                //  voxel(NULL),
								  dens_(nullptr),
                                  vX_(nullptr),
                                  vY_(nullptr),
                                  vZ_(nullptr)
{
}

VolumeInstance::~VolumeInstance()
{
	delete[] dens_;
	dens_ = nullptr;
	delete[] vX_;
	vX_ = nullptr;
	delete[] vY_;
	vY_ = nullptr;
	delete[] vZ_;
	vZ_ = nullptr;
	delete[] tempK_;
	delete[] impl_surf_;
	delete[] press_;
	delete[] f_conf;
	delete[] Y_elaps;
	//delete[] voxel;
}

void VolumeInstance::allocate_memory(unsigned long long int size)
{
//	voxel = new Voxel[size];
	dens_ = new double[size];
	vX_ = new double[size*2];
	vY_ = new double[size*2];
	vZ_ = new double[size*2];
	tempK_ = new double[size];
	press_ = new double[size];
	impl_surf_ = new int[size];
	f_conf = new double[size];
	Y_elaps = new double[size];


	visc_ = AIR_VISC;
	diff_ = AIR_DIFF;
	dt_ = DEFAULT_time_step;
}

void VolumeInstance::add_source(double * x, double * a, unsigned long int n)
{
	unsigned long int i, j, k;

	for (i = 0; i <= n + 1; i++)
	{
		for (j = 0; j <= n + 1; j++)
		{
			for (k = 0; k <= n + 1; k++)
			{
				x[IDX(i, j, k, n)] += a[IDX(i, j, k, n)];
			}
		}
	}
}

void VolumeInstance::add_constant(double* x, double a, unsigned long int n)
{
	unsigned long int i, j, k;

	for (i = 0; i <= n + 1; i++)
	{
		for (j = 0; j <= n + 1; j++)
		{
			for (k = 0; k <= n + 1; k++)
			{
				x[IDX(i, j, k, n)] += a;
			}
		}
	}
}

void VolumeInstance::draw_sphere(unsigned long int n)
{
	long long int i, j, k;
	const unsigned long int half_n = static_cast<int>(n / 2);
	max_ = 0; min_ = 0;

	for (i = 0; i <= n + 1; i++)
	{
		for (j = 0; j <= n + 1; j++)
		{
			for (k = 0; k <= n + 1; k++)
			{
				const double density_value = ((i <= half_n) ? i : (n + 2 - i)) * ((j <= half_n) ? j : (n + 2 - j)) *(
					(k <= half_n) ? k : (n + 2 - k)); 
				dens_[IDX(i, j, k, n)] = density_value;
				vX_[IDX(i, j, k, n)] = 0;
				vY_[IDX(i, j, k, n)] = DEFAULT_high_velocity;
				vZ_[IDX(i, j, k, n)] = DEFAULT_avg_velocity;
				tempK_[IDX(i, j, k, n)] = KELVINIZE(T_AIR);
				f_conf[IDX(i, j, k, n)] = 0;
				press_[IDX(i, j, k, n)] = 0;
				Y_elaps[IDX(i, j, k, n)] = 0;
				/*
				if(density_value > max_)
				{
					max_ = density_value;
				}
				if (density_value < min_)
				{
					min_ = density_value;
				}
				*/
			}
		}
	}


	add_fuel(n);
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
				dens_[IDX(i, j, k, n)] = 0;
				vX_[IDX(i, j, k, n)] = 0;
				vY_[IDX(i, j, k, n)] = 0;
				vZ_[IDX(i, j, k, n)] = 0;
				tempK_[IDX(i, j, k, n)] = KELVINIZE(T_AIR);
				f_conf[IDX(i, j, k, n)] = 0;
				press_[IDX(i, j, k, n)] = 0;
				Y_elaps[IDX(i, j, k, n)] = 0;
			}
		}
	}
}

void VolumeInstance::copy(unsigned long int n, VolumeInstance* from_volume)
{
	
	unsigned long int i, j, k;
	max_ = from_volume->max_; min_ = from_volume->min_;
	visc_ = from_volume->visc_; diff_ = from_volume->diff_; dt_ = from_volume->dt_;

	for (i = 0; i <= n + 1; i++)
	{
		for (j = 0; j <= n + 1; j++)
		{
			for (k = 0; k <= n + 1; k++)
			{
				dens_[IDX(i, j, k, n)] = from_volume->dens_[IDX(i, j, k, n)];
				vX_[IDX(i, j, k, n)] = from_volume->vX_[IDX(i, j, k, n)];
				vY_[IDX(i, j, k, n)] = from_volume->vY_[IDX(i, j, k, n)];
				vZ_[IDX(i, j, k, n)] = from_volume->vZ_[IDX(i, j, k, n)];
				press_[IDX(i, j, k, n)] = from_volume->press_[IDX(i, j, k, n)];
				tempK_[IDX(i, j, k, n)] = from_volume->tempK_[IDX(i, j, k, n)];
				impl_surf_[IDX(i, j, k, n)] = from_volume->impl_surf_[IDX(i, j, k, n)];
				f_conf[IDX(i, j, k, n)] = from_volume->f_conf[IDX(i, j, k, n)];

				Y_elaps[IDX(i, j, k, n)] = from_volume->Y_elaps[IDX(i, j, k, n)];

			}
		}
	}
}

void VolumeInstance::add_fuel(unsigned long int n)
{
	long int half_n = static_cast<long int>(n / 2);
	long long int i, j, k;

	for (i = 0; i <= n+1; i++)
	{
		for (j = 0; j <= n+1; j++)
		{
			for (k = 0; k <= n+1; k++)
			{
				if (STD_1g_BLUE_CORE(VXL_TO_M(i- half_n), VXL_TO_M(j- half_n), VXL_TO_M(k)- WICK_H * 2) <= 1.0) {
					Y_elaps[IDX(i,j,k, n)] = 1;
				//	cout << "++++++++++++++++++" << endl;
				}

			//	cout << " " << VXL_TO_M(i-half_n ) << " " << VXL_TO_M(j-half_n ) << " " << VXL_TO_M(k) + 0.01 << "     " << STD_1g_BLUE_CORE(VXL_TO_M(i-half_n), VXL_TO_M(j-half_n), VXL_TO_M(k)+0.01)  << endl;
			}
		}
	}
}
