#pragma once

enum voxelType { boundary, wick, air }; //, gasFuel, gasProducts};

class Voxel
{
private:
	
public:
	double vX, vY, vZ, // velocities
		tempC, // temperature in C
		dens, // density
		p = 9.8, // pressure
		Rspec, // specific gas constant
		C25H52, O2, // agents of reaction
		C, H2O, CO, CO2; // products of reaction

	// density of air can be calculated using the ideal gas law: 
	// dens = p / (Rspec*T)
	//Rspec - specific gas constant for dry air

	Voxel();
	~Voxel();
	double get_dens();
	double get_abs_temp();
};

