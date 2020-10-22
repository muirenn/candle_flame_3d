#include "FullAnimation.h"

void FullAnimation::init() {
    iter = 0;
    write_init_data_to_file();

    frame_size_ = (n_ + 2) * (n_ + 2) * (n_ + 2); // +2 goes for borders at each side

    frame.allocate_memory(frame_size_);
    frame_prev.allocate_memory(frame_size_);

    unsigned long long int d_length = time_ * frame_size_;
    qty_to_display = new(nothrow) double[d_length];
    if (qty_to_display == nullptr) {
        cerr << "MEMORY ALLOCATION ERROR" << endl;
    }

    frame.hardcode_init_values(n_);

    frame_prev.fill_with_zeros(n_);

    minTotal_ = 0;
    maxTotal_ = 0;
}

void FullAnimation::run() {
    // define/track the blue core
    frame.add_fuel(n_);

    unsigned int t;

    for (t = 0; t < time_; t++) {
        // fbuoy = a (T −Tair)(0,0,1)
        apply_buoyancy_forces();
        // apply_confinement_forces();

        /* Velocity step: ∂u/∂t=-(u∙∇)u+ν∇^2 u+f */

        // ν∇^2 u
        diffuse(frame.dt_, frame.visc_, velocity_x, frame_prev.vX_, frame.vX_, lin_itr_, n_);
        diffuse(frame.dt_, frame.visc_, velocity_y, frame_prev.vY_, frame.vY_, lin_itr_, n_);
        diffuse(frame.dt_, frame.visc_, velocity_z, frame_prev.vZ_, frame.vZ_, lin_itr_, n_);

        project(frame_prev.vX_, frame_prev.vY_, frame_prev.vZ_, frame.vX_, frame.vY_, lin_itr_, n_);

        // -(u∙∇)u for each part of u
        advect(frame.dt_, velocity_x, frame.vX_, frame_prev.vX_, frame_prev.vX_, frame_prev.vY_, frame_prev.vZ_, n_);
        advect(frame.dt_, velocity_y, frame.vY_, frame_prev.vY_, frame_prev.vX_, frame_prev.vY_, frame_prev.vZ_, n_);
        advect(frame.dt_, velocity_z, frame.vZ_, frame_prev.vZ_, frame_prev.vX_, frame_prev.vY_, frame_prev.vZ_, n_);

        project(frame.vX_, frame.vY_, frame.vZ_, frame_prev.vX_, frame_prev.vY_, lin_itr_, n_);

        // TODO add pressure (correction) term:
        // u = u*-∆t∇p/ρ
        // solve ∇∙(∇p/ρ)=∇∙u^*/∆t to find the pressure needed to update the above equation

        /* Density step: ∂ρ/∂t=-(u∙∇)ρ+κ∇^2 ρ+S*/
       diffuse(frame.dt_, frame.diff_, density, frame_prev.dens_, frame.dens_, lin_itr_, n_);
       advect(frame.dt_, density, frame.dens_, frame_prev.dens_, frame.vX_, frame.vY_, frame.vZ_, n_);

        write_results();

        iter++;
        frame.add_fuel(n_);
        frame_prev.copy(n_, &frame);
    }
}

void FullAnimation::write_init_data_to_file() const {
    ofstream blender_file, log_file;
    blender_file.open(filename_, ios::out | ios::binary | ios::trunc);
    log_file.open(log_filename_, ios::out | ios::trunc);
    blender_file.write(reinterpret_cast<const char *>(&n_), sizeof(n_));
    blender_file.write(reinterpret_cast<const char *>(&n_), sizeof(n_));
    blender_file.write(reinterpret_cast<const char *>(&n_), sizeof(n_));
    blender_file.write(reinterpret_cast<const char *>(&time_), sizeof(time_));
    log_file << "Sides: " << n_ << " units;" <<  endl << "Time: " << time_ << " steps." << endl;
    blender_file.close();
    log_file.close();
}

void FullAnimation::write_animation_to_file() {
    ofstream blender_file, log_file;
    blender_file.open(filename_, ios::out | ios::binary | ios::app);
    log_file.open(log_filename_, ios::out | ios::app);

    unsigned long int i, j, k, l;

    for (l = 0; l < time_; l++) {
        log_file << endl << "===== TIMESTAMP " << l << " =====" << endl;
        for (k = 1; k <= n_; k++) {
            for (j = 1; j <= n_; j++) {
                for (i = 1; i <= n_; i++) {
                    double qty_unnormalized = qty_to_display[l * frame_size_ + IDX(i, j, k, n_)];
                    float qty_normal = static_cast<float>((qty_unnormalized - minTotal_) / (maxTotal_ - minTotal_));
                    blender_file.write(reinterpret_cast<const char *>(&qty_normal), sizeof(qty_normal));
                    log_file << std::setw(10) << qty_normal << " ";
                }
                log_file << endl;
            }
            log_file << "-----------------------------------------" << endl;
        }
    }


    blender_file.close();
    log_file.close();
}

// The temperature affects the fluid velocity as hot gases tend to rise due to buoyancy
void FullAnimation::apply_buoyancy_forces() {
    unsigned long int i, j, k;
    for (i = 0; i <= n_ + 1; i++) {
        for (j = 0; j <= n_ + 1; j++) {
            for (k = 0; k <= n_ + 1; k++) {
                frame.vZ_[IDX(i, j, k, n_)] += ALPHA * frame.dens_[IDX(i, j, k, n_)];
            }
        }
    }
}

void FullAnimation::apply_confinement_forces() {
    unsigned long int i, j, k;
    double dx, dy, dz, len;
    double vorticity = 10;

    for (i = 2; i <= n_ - 2; i++) {
        for (j = 2; j <= n_ - 2; j++) {
            for (k = 2; k <= n_ - 2; k++) {
                double absCurlX =  abs(curl(frame.vX_,frame.vY_,frame.vZ_,i - 1,j,k,n_)) - abs(curl(frame.vX_,frame.vY_,frame.vZ_,i + 1,j,k,n_));
                double absCurlY =  abs(curl(frame.vX_,frame.vY_,frame.vZ_,i,j - 1,k,n_)) - abs(curl(frame.vX_,frame.vY_,frame.vZ_,i,j + 1,k,n_));
                double absCurlZ = abs(curl(frame.vX_,frame.vY_,frame.vZ_,i,j,k - 1,n_)) - abs(curl(frame.vX_,frame.vY_,frame.vZ_,i,j,k + 1,n_));

                // directions to X, Y, Z
                dx = absCurlY + absCurlZ;
                dy = absCurlX + absCurlZ;
                dz = absCurlX + absCurlY;

                len = sqrt(dx*dx + dy*dy + dz*dz) + 1e-5f;
                dx = vorticity/len*dx;
                dy = vorticity/len*dy;
                dz = vorticity/len*dz;

                double curlXYZ = curl(frame.vX_,frame.vY_,frame.vZ_,i,j,k,n_);

                frame.vX_[IDX(i,j,k,n_)] = frame.vX_[IDX(i,j,k,n_)] + frame.dt_*curlXYZ*dx;
                frame.vY_[IDX(i,j,k,n_)] = frame.vY_[IDX(i,j,k,n_)] + frame.dt_*curlXYZ*dy;
                frame.vZ_[IDX(i,j,k,n_)] = frame.vZ_[IDX(i,j,k,n_)] + frame.dt_*curlXYZ*dz;
            }
        }
    }


}

void FullAnimation::write_results() {
    unsigned long int i, j, k;


    for (i = 1; i <= n_; i++) {
        for (j = 1; j <= n_; j++) {
            for (k = 1; k <= n_; k++) {
                double res = frame.dens_[IDX(i, j, k, n_)];
                qty_to_display[iter * frame_size_ + IDX(i, j, k, n_)] = res; // TODO check order
                if (res < minTotal_) {
                    minTotal_ = res;
                }
                if (res > maxTotal_) {
                    maxTotal_ = res;
                }
            }
        }
    }
}

FullAnimation::FullAnimation() : qty_to_display(nullptr), time_(N_FRAMES),
                                 n_(M_TO_VXL(N_SIDE)),
                                 lin_itr_(LIN_ITERS),
                                 frame_size_(0), iter(0),
                                 minTotal_(0),
                                 maxTotal_(0) {
}

FullAnimation::FullAnimation(const unsigned int time, const unsigned int n, const string filename, const string log_filename)
        : FullAnimation() {
    time_ = time;
    n_ = n;
    filename_ = filename;
    log_filename_ = log_filename;
}


FullAnimation::~FullAnimation() {
    delete[] qty_to_display;
}
