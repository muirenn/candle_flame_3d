#pragma once

/* User-defined & user-adjustable constants */
#define TIME 10       // [s]		real physical time of the animation
#define N_SIDE 0.04   // [m]		real physical cube side size
#define SPACING 0.001 // [m/vxl]  uniform spacing
#define TIME_STEP 0.1 // [s/frame]
#define WICK_D 0.001  // [m]		wick diameter
#define WICK_H 0.005  // [m]		wick height
#define LIN_ITERS 20  //			number of iterations for linear solver

/* Known physical quantities & assumed constants */
#define RO_F 1         // [kg/m^3]	density of the gaseous fuel - gaseous paraffin
#define RO_H 0.01      // [kg/m^3]	density f hot gasesous products
#define RO_S 770       // [kg/m^3]	density of the solid fuel (melted paraffin)
#define EPS_F 16       //  vorticity confinement parameter for gaseous fuel
#define EPS_H 60       //	vorticity confinement parameter for hot gaseous products
#define SOURCE 0.1     // [m/s]
#define T_IGN 600      // [degr C]	ignition temperature
#define T_MAX 1400     // [degr C]	maximum flame temperature
#define T_PRF_MELT 55  // [degr C]	melting point of paraffin
#define T_PRF_BOIL 400 // [degr C]	boiling point of paraffin
#define T_AIR 21       // [degr C] temerature of the surrounding air
#define FL_SPEED 0.5   // [m/s]	or 0.161: flame reaction speed - property of the fuel
#define ALPHA 10       // [m/Ks^2]	positive constant
#define CONST_K 1      //			positive constant
#define C_T 3000       // [K/s]	cooling constant
#define CONST_S 0.5    // [m/s]    property of the fuel, TODO adjust
#define V_S 0.000029   // [kg/s]

#define AIR_VISC_T_ROOM 0.000015  // [m^2/s]	air viscosity
#define CO2_VISC_T_IGN 0.000062   //
#define CO2_VISC_T_ROOM 0.000008  //
#define AIR_DIFF_T_ROOM 0.00002   // [m^2/s]  air diffusivity	0.00002	L^2*T^(-1); L = initial size of gas blob or field
#define CO2_DIFF_T_IGN 0.000085   //
#define CO2_DIFF_T_ROOM 0.000011  //
#define STD_VORTICITY 10          //
#define DEFAULT_avg_density 5     //
#define DEFAULT_high_density 20   //
#define DEFAULT_avg_velocity 0.8  // [m/s?]
#define DEFAULT_high_velocity 1.2 //

/* Molecular weights of chemical compounds [g/mol] */
#define MW_C25H52 352.68 // paraffin wax (pentacosane)
#define MW_O2 32         // oxygen
#define MW_CO2 44.01     // carbon dioxide
#define MW_H2O 18.02     // water
#define MW_CO 28.01      // carbon monoxide
#define MW_C 12.01       // carbon

/* Parameters of the fuel region in the flames. In this case: 1g (normal gravity), measured & approximated in (SRC TODO) */
#define CORE_H_R_1g 0.003        // [m] flame core horizontal radius for a candle in normal gravity
#define CORE_V_R_1g 0.005        // [m] flame core vertical radius for a candle in normal gravity
#define DARK_REGION_H_R_1g 0.002 // [m] dark region (no-oxygen) horizontal radius in normal gravity
#define DARK_REGION_V_R_1g 0.005 // [m] dark region (no-oxygen) vertical radius in normal gravity
#define DARK_REGION_SHIFT 0.0    // [m]

/* Conversion of units */
#define KELVINIZE(t) ((t) + 273)
#define UNKELVINIZE(t) ((t)-273)
#define CUBE_SIDE (N_SIDE / SPACING)
#define VXL_TO_M(x) ((x) * (SPACING))
#define M_TO_VXL(x) ((x) / SPACING)
#define V_F (V_S + (RO_S / RO_F - 1) * CONST_S)
#define N_FRAMES (TIME / TIME_STEP)
#define S_TO_FRAME(x) ((x) / TIME_STEP)

#define M2perS_TO_VXL2perFRAME(x) ((M_TO_VXL(x)) * (M_TO_VXL(x)) / (S_TO_FRAME(x)))

/* ID of element (i,j,k) in the virtual cube (3-dimensional array) with side N */
#define IDX(i, j, k, N) ((k) + (j) * ((N) + 2) + (i) * ((N) + 2) * ((N) + 2)) // i === z for loop with i->j->k

/* Checks whether the point with coordinates (i,j,k) lies within the blue core of the flame given normal gravity
 * If result is <=1 then the point lies within the ellipsoid */
#define STD_1g_BLUE_CORE(i, j, k)                            \
    (((i) * (i) + (j) * (j)) / (CORE_H_R_1g * CORE_H_R_1g) + \
     (k) * (k) / (CORE_V_R_1g * CORE_V_R_1g))

/* Checks whether the point with coordinates (i,j,k) lies within the dark (non-luminous) zone of the flame given normal
 * gravity
 * If result is <=1 then the point lies within the ellipsoid */
#define STD_1g_DARK_ZONE(i, j, k)                                          \
    (((i) * (i) + (j) * (j)) / (DARK_REGION_H_R_1g * DARK_REGION_H_R_1g) + \
     (k + DARK_REGION_SHIFT) * (k + DARK_REGION_SHIFT) / (DARK_REGION_V_R_1g * DARK_REGION_V_R_1g))

constexpr const char *bvox_filename_DEFAULT = "candleflame.bvox";
#define LOG_NAME_LENGTH 32
