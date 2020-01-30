#include "stdafx.h"
#include "Voxel.h"


Voxel::Voxel()
{
}


Voxel::~Voxel()
{
}

double Voxel::get_dens()
{
	return p / (Rspec * get_abs_temp());
}

double Voxel::get_abs_temp()
{
	return tempC + 273.0;
}
