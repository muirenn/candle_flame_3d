#pragma once


/* Known physical quantities & assumed constants */
#define TIME 5			// [s]		real physical time of the animation
#define N_SIDE 0.04		// [m]		cube side size
#define RO_F 1			// [kg/m^3]	density of the gaseous fuel - gaseous paraffin
#define RO_H 0.01		// [kg/m^3]	density f hot gasesous products
#define RO_S 770		// [kg/m^3]	density of the solid fuel (melted paraffin)
#define EPS_F 16		//  vorticity confinement parameter for gaseous fuel
#define EPS_H 60		//	vorticity confinement parameter for hot gaseous products
#define SOURCE 0.1		// [m/s]
#define SPACING 0.0005	// [m/vxl]  uniform spacing
#define T_IGN 600		// [degr C]	ignition temperature
#define T_MAX 1400		// [degr C]	maximum flame temperature
#define T_PRF_MELT 55	// [degr C]	melting point of paraffin
#define T_PRF_BOIL 400	// [degr C]	boiling point of paraffin
#define T_AIR 21		// [degr C] temerature of the surrounding air
#define FL_SPEED 0.5	// [m/s]	or 0.161: flame reaction speed - property of the fuel	
#define ALPHA 3		// [m/Ks^2]	positive constant
#define TIME_STEP 0.1	// [s/frame]
#define CONST_K 1		//			positive constant
#define C_T 3000		// [K/s]	cooling constant
#define AIR_VISC_T_ROOM 0.000015// [m^2/s]	air viscosity
#define CO2_VISC_T_IGN 0.000062
#define CO2_VISC_T_ROOM 0.000008
#define AIR_DIFF_T_ROOM 0.00002      // [m^2/s]  air diffusivity	0.00002	L^2*T^(-1); L = initial size of gas blob or field
#define CO2_DIFF_T_IGN 0.000085
#define CO2_DIFF_T_ROOM 0.000011

#define CONST_S 0.5		// [m/s]    property of the fuel, TODO adjust
#define V_S 0.000029	// [kg/s]

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
#define STD_1g_DARK_ZONE(i,j,k) (((i)*(i)+(j)*(j))/(0.0025*0.00025)+(k)*(k)/(0.004*0.004))
#define V_F (V_S+(RO_S/RO_F - 1)*CONST_S)
#define N_FRAMES (TIME/TIME_STEP)
#define S_TO_FRAME(x) ((x)/TIME_STEP)


#define KELVINIZE(t) ((t) + 273)
#define UNKELVINIZE(t) ((t) - 273)
#define CUBE_SIDE (N_SIDE/SPACING)
#define VXL_TO_M(x) ((x)*(SPACING))
#define M_TO_VXL(x) ((x)/SPACING)
#define M2perS_TO_VXL2perFRAME(x) ((M_TO_VXL(x))*(M_TO_VXL(x))/(S_TO_FRAME(x)))


#define IDX(i,j,k,N) ((k)+(j)*((N)+2)+(i)*((N)+2)*((N)+2)) // i === z for loop with i->j->k


// to rename
#define DEFAULT_cube_side_size 50

#define DEFAULT_avg_density 5
#define DEFAULT_high_density 200
#define DEFAULT_avg_velocity 0.8// [m/s?]
#define DEFAULT_high_velocity 1.2
constexpr const char* bvox_filename_DEFAULT = "candleflame.bvox";
