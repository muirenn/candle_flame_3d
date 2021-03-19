#include "VolumeInstance.h"
#include "Math_Solver.h"

VolumeInstance::VolumeInstance() : visc_(AIR_VISC_T_ROOM),
                                   diff_(AIR_DIFF_T_ROOM), dt_(TIME_STEP), min_(0), max_(0),
                                   dens(nullptr),
                                   vX(nullptr),
                                   vY(nullptr),
                                   vZ(nullptr) {
}

VolumeInstance::~VolumeInstance() {
    delete[] dens;
    dens = nullptr;
    delete[] vX;
    vX = nullptr;
    delete[] vY;
    vY = nullptr;
    delete[] vZ;
    vZ = nullptr;
}

void VolumeInstance::allocate_memory(unsigned long long int size) {
    dens = new double[size];
    vX = new double[size * 2];
    vY = new double[size * 2];
    vZ = new double[size * 2];

    visc_ = CO2_VISC_T_IGN;
    diff_ = AIR_DIFF_T_ROOM;
    dt_ = TIME_STEP;
}

void VolumeInstance::add_source(double *x, double *a, unsigned long int n) {
    unsigned long int i, j, k;

    for (i = 0; i <= n + 1; i++) {
        for (j = 0; j <= n + 1; j++) {
            for (k = 0; k <= n + 1; k++) {
                x[IDX(i, j, k, n)] += a[IDX(i, j, k, n)];
            }
        }
    }
}

void VolumeInstance::add_constant(double *x, double a, unsigned long int n) {
    unsigned long int i, j, k;

    for (i = 0; i <= n + 1; i++) {
        for (j = 0; j <= n + 1; j++) {
            for (k = 0; k <= n + 1; k++) {
                x[IDX(i, j, k, n)] += a;
            }
        }
    }
}

void VolumeInstance::fill_with_zeros(unsigned long int n) {
    unsigned long int i, j, k;
    max_ = 0;
    min_ = KELVINIZE(T_MAX);

    for (i = 0; i <= n + 1; i++) {
        for (j = 0; j <= n + 1; j++) {
            for (k = 0; k <= n + 1; k++) {
                dens[IDX(i, j, k, n)] = 0;
                vX[IDX(i, j, k, n)] = 0;
                vY[IDX(i, j, k, n)] = 0;
                vZ[IDX(i, j, k, n)] = 0;
            }
        }
    }
}

void VolumeInstance::hardcode_init_values(unsigned long int n) {
    unsigned long int i, j, k;
    max_ = 0;
    min_ = KELVINIZE(T_MAX);
    long int half_n = static_cast<long int>(n / 2);

    for (i = 0; i <= n + 1; i++) {
        for (j = 0; j <= n + 1; j++) {
            for (k = 0; k <= n + 1; k++) {
                dens[IDX(i, j, k, n)] = 0;
                vX[IDX(i, j, k, n)] = 0;
                vY[IDX(i, j, k, n)] = 0;
                vZ[IDX(i, j, k, n)] = 0;

                // If the current voxel is within the area of estimated "blue core", we assume it has the temperature of ignition.
                // Otherwise, it has air temperature
                long i0 = (i - half_n), j0 = (j - half_n);
                if (fuel_zone(i0,j0,k)) {
                    dens[IDX(i, k, j, n)] = DEFAULT_high_density;
                } else {
                    dens[IDX(i, k, j, n)] = 0;
                }
            }
        }
    }
}

void VolumeInstance::copy(unsigned long int n, VolumeInstance *from_volume) {

    unsigned long int i, j, k;
    max_ = from_volume->max_;
    min_ = from_volume->min_;
    visc_ = from_volume->visc_;
    diff_ = from_volume->diff_;
    dt_ = from_volume->dt_;

    for (i = 0; i <= n + 1; i++) {
        for (j = 0; j <= n + 1; j++) {
            for (k = 0; k <= n + 1; k++) {
                dens[IDX(i, j, k, n)] = from_volume->dens[IDX(i, j, k, n)];
                vX[IDX(i, j, k, n)] = from_volume->vX[IDX(i, j, k, n)];
                vY[IDX(i, j, k, n)] = from_volume->vY[IDX(i, j, k, n)];
                vZ[IDX(i, j, k, n)] = from_volume->vZ[IDX(i, j, k, n)];
            }
        }
    }
}

void VolumeInstance::add_fuel(unsigned long int n) {
    long int half_n = static_cast<long int>(n / 2);
    long long int i, j, k;

    for (i = 0; i <= n + 1; i++) {
        for (j = 0; j <= n + 1; j++) {
            for (k = 0; k <= n + 1; k++) {
                long i0 = (i - half_n), j0 = (j - half_n);
                if (fuel_zone(i0,j0,k)) {

                    dens[IDX(i, j, k, n)] += DEFAULT_high_density;
                }
            }
        }
    }
}
