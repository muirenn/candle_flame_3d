#pragma once


/* Known physical quantities & assumed constants */
#define N_SIDE 0.05		// [m]		cube side size
#define RO_F 1			// [kg/m^3]	density of the gaseous fuel - gaseous paraffin
#define RO_H 0.01		// [kg/m^3]	density f hot gasesous products
#define RO_S 770		// [kg/m^3]	density of the solid fuel (melted paraffin)
#define EPS_F 16		//			vorticity confinement parameter for gaseous fuel
#define EPS_H 60		//			vorticity confinement parameter for hot gaseous products
#define SOURCE 0.1		// [m/s]
#define SPACING 0.001	// [m/voxel]uniform spacing
#define T_IGN 260		// [degr C]	ignition temperature TODO convert to K
#define T_MAX 1400		// [degr C]	maximum flame temperature
#define T_PRF_MELT 55	// [degr C]	melting point of paraffin
#define T_PRF_BOIL 400	// [degr C]	boiling point of paraffin
#define T_AIR 21		// [degr C] temerature of the surrounding air
#define FL_SPEED 0.5	// [m/s]	or 0.161: flame reaction speed - property of the fuel	
#define ALPHA 0.15		// [m/Ks^2]	positive constant
#define TIME_STEP 0.2	// [s]	
#define CONST_K 1		//			positive constant
#define C_T 3000		// [K/s]	cooling constant
#define AIR_VISC 0.000015// [m^2/s]	air viscosity
#define AIR_DIFF 0.00002// [m^2/s] air diffusivity		L^2*T^(-1); L = initial size of gas blob or field



#define WICK_D 0.001	// [m]		wick diameter
#define WICK_H 0.005	// [m]		wick height
#define LIN_ITERS 20	//			number of iterations for linear solver

/* Molecular weights of chemical compounds [g/mol] */
#define MW_C25H52 352.68	// paraffin wax (pentacosane)
#define MW_O2 32			// oxygen	
#define MW_CO2 44.01		// carbon dioxide
#define MW_H2O 18.02		// water
#define MW_CO 28.01			// carbon monoxide
#define MW_C 12.01			// carbon


#define STD_1g_WICK_LENGTH 0.004 // [m]
#define STD_1g_BLUE_CORE(i,j,k) (((i)*(i)+(j)*(j))/(0.003*0.003)+(k)*(k)/(0.005*0.005))


#define KELVINIZE(t) ((t) + 273)
#define UNKELVINIZE(t) ((t) - 273)
#define CUBE_SIDE (N_SIDE/SPACING)
#define VXL_TO_M(x) ((x)*(SPACING))
#define M_TO_VXL(x) (x/SPACING)




#define IDX(i,j,k,N) ((k)+(j)*((N)+2)+(i)*((N)+2)*((N)+2)) // i === z for loop with i->j->k

#define CANDLE_TO_BVOX(i,j,k,n) () // translate coordinates to candle coordinate system

#define TO_CANDLE(i,j,k,nd2,n) (IDX(i-nd2,j-nd2,k,n))



// to rename
#define DEFAULT_cube_side_size 50
#define DEFAULT_time_step 0.2 // T

#define DEFAULT_avg_density 5
#define DEFAULT_high_density 20
#define DEFAULT_avg_velocity 0.2// [m/s?]
#define DEFAULT_high_velocity 0
constexpr const char* bvox_filename_DEFAULT = "candleflame.bvox";
