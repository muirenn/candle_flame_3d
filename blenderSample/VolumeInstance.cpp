#include "VolumeInstance.h"

VolumeInstance::VolumeInstance() : visc_(AIR_VISC_T_ROOM),
                                   diff_(AIR_DIFF_T_ROOM), dt_(TIME_STEP), min_(0), max_(0),
                                   dens_(nullptr),
                                   vX_(nullptr),
                                   vY_(nullptr),
                                   vZ_(nullptr) {
}

VolumeInstance::~VolumeInstance() {
    delete[] dens_;
    dens_ = nullptr;
    delete[] vX_;
    vX_ = nullptr;
    delete[] vY_;
    vY_ = nullptr;
    delete[] vZ_;
    vZ_ = nullptr;
}

void VolumeInstance::allocate_memory(unsigned long long int size) {
    dens_ = new double[size];
    vX_ = new double[size * 2];
    vY_ = new double[size * 2];
    vZ_ = new double[size * 2];

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
                dens_[IDX(i, j, k, n)] = 0;
                vX_[IDX(i, j, k, n)] = 0;
                vY_[IDX(i, j, k, n)] = 0;
                vZ_[IDX(i, j, k, n)] = 0;
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
                dens_[IDX(i, j, k, n)] = 0;
                vX_[IDX(i, j, k, n)] = 0;
                vY_[IDX(i, j, k, n)] = 0;
                vZ_[IDX(i, j, k, n)] = 0;

                // If the current voxel is within the area of estimated "blue core", we assume it has the temperature of ignition.
                // Otherwise, it has air temperature
                // For simplification, STG_1g_BLUE_CORE is a sphere in the center of the scene
                long i0 = labs(i - half_n), j0 = labs(j - half_n);
                if ((STD_1g_BLUE_CORE(VXL_TO_M(i0), VXL_TO_M(j0), VXL_TO_M(k) - WICK_H * 2) <= 1.0) &&
                    (STD_1g_DARK_ZONE(VXL_TO_M(i0), VXL_TO_M(j0), VXL_TO_M(k) - WICK_H * 2) > 1.0)) {
                    dens_[IDX(i, j, k, n)] = DEFAULT_high_density;
                } else {
                    dens_[IDX(i, j, k, n)] = 0;
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
                dens_[IDX(i, j, k, n)] = from_volume->dens_[IDX(i, j, k, n)];
                vX_[IDX(i, j, k, n)] = from_volume->vX_[IDX(i, j, k, n)];
                vY_[IDX(i, j, k, n)] = from_volume->vY_[IDX(i, j, k, n)];
                vZ_[IDX(i, j, k, n)] = from_volume->vZ_[IDX(i, j, k, n)];
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
                long i0 = labs(i - half_n), j0 = labs(j - half_n);
                if ((STD_1g_BLUE_CORE(VXL_TO_M(i0), VXL_TO_M(j0), VXL_TO_M(k) - WICK_H * 2) <= 1.0) &&
                    (STD_1g_DARK_ZONE(VXL_TO_M(i0), VXL_TO_M(j0), VXL_TO_M(k) - WICK_H * 2) > 1.0)) {

                    dens_[IDX(i, j, k, n)] += DEFAULT_high_density;
                }
            }
        }
    }
}
