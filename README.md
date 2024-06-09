# Modelling and Visualisation of a Candle Flame

**Abstract**. This thesis delves into the creation of a 3D simulation of a candle flame, drawing inspiration from Jos Stam’s groundbreaking works in fluid dynamics. Stam’s algorithms, initially developed for 2D simulations, are extended to a 3D domain, enabling the detailed visualization of candle flames. The simulation employs computational fluid dynamics (CFD) principles, capturing the intricate flow patterns and temperature distribution within the flame. The governing equations are solved using the Gauss-Seidel relaxation method, implemented in the C++ programming language. The resulting simulation is visualized in Blender 2.79, showcasing the intricate dynamics of the candle flame. This work demonstrates the potential of CFD in 3D simulations and paves the way for further exploration of fluid dynamics phenomena in a 3D environment.

**Keywords**: Computational Fluid Dynamics, Gauss-Seidel Relaxation, Volumetric Datacubes


## Features

- **3D Fluid Simulation**: Implements the Navier-Stokes equations to simulate fluid flow in three dimensions.
- **Density and Velocity Fields**: Updates density and velocity values at each time step based on environmental factors and external forces.
- **Volumetric Rendering**: Uses a volumetric ray tracer for rendering the fluid in Blender.
- **Gauss-Seidel Relaxation**: Employs an iterative method to solve the linear systems for diffusion.
- **Linear Interpolation**: Uses linear interpolation for advection steps to update fluid properties.
- **Boundary Handling**: Manages the boundaries of the simulation domain to ensure realistic fluid behavior.

## Algorithm

The fluid solver operates in a loop, updating the velocity and density fields at each time step. The core steps in the solver include:

1. **Initialization**: Set up initial density and velocity values.
2. **Velocity Step**: 
    - Apply buoyancy and confinement forces.
    - Diffuse the velocity fields using Gauss-Seidel relaxation.
    - Project the velocity field to ensure it is divergence-free.
    - Advect the velocity fields.
3. **Density Step**:
    - Diffuse the density field.
    - Advect the density field using the updated velocity field.
4. **Boundary Conditions**: Apply boundary conditions to maintain fluid properties within the simulation domain.

## Implementation

### Diffusion Step

The diffusion step uses Gauss-Seidel relaxation to iteratively solve the system of linear equations representing the diffusion of velocity and density fields. This method ensures fast convergence for large, sparse systems.

```cpp
void diffuse(int N, int b, double *x, double *x0, double diff, double dt) {
    double a = dt * diff * N * N;
    for (int k = 0; k < 20; k++) {
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i-1, j)] + x[IX(i+1, j)] + x[IX(i, j-1)] + x[IX(i, j+1)])) / (1 + 4 * a);
            }
        }
        set_bnd(N, b, x);
    }
}
```

### Advection Step

The advection step updates the properties of the fluid by moving them along the velocity field. Linear interpolation is used to compute the new values.

```cpp
void advect(int N, int b, double *d, double *d0, double *u, double *v, double dt) {
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            double x = i - dt * u[IX(i, j)];
            double y = j - dt * v[IX(i, j)];
            if (x < 0.5) x = 0.5; if (x > N + 0.5) x = N + 0.5;
            if (y < 0.5) y = 0.5; if (y > N + 0.5) y = N + 0.5;
            int i0 = (int) x, i1 = i0 + 1;
            int j0 = (int) y, j1 = j0 + 1;
            double s1 = x - i0, s0 = 1 - s1;
            double t1 = y - j0, t0 = 1 - t1;
            d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) + s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
        }
    }
    set_bnd(N, b, d);
}
```

### Projection Step
The projection step ensures that the velocity field is divergence-free, which is essential for maintaining the incompressibility of the fluid.

```cpp
void project(int N, double *u, double *v, double *p, double *div) {
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            div[IX(i, j)] = -0.5 * (u[IX(i+1, j)] - u[IX(i-1, j)] + v[IX(i, j+1)] - v[IX(i, j-1)]) / N;
            p[IX(i, j)] = 0;
        }
    }
    set_bnd(N, 0, div);
    set_bnd(N, 0, p);

    for (int k = 0; k < 20; k++) {
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                p[IX(i, j)] = (div[IX(i, j)] + p[IX(i-1, j)] + p[IX(i+1, j)] + p[IX(i, j-1)] + p[IX(i, j+1)]) / 4;
            }
        }
        set_bnd(N, 0, p);
    }

    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            u[IX(i, j)] -= 0.5 * N * (p[IX(i+1, j)] - p[IX(i-1, j)]);
            v[IX(i, j)] -= 0.5 * N * (p[IX(i, j+1)] - p[IX(i, j-1)]);
        }
    }
    set_bnd(N, 1, u);
    set_bnd(N, 2, v);
}
```

### Visualization
The fluid simulation is visualized using a volumetric ray tracer in Blender. The solver outputs density and velocity fields, which are then rendered to produce realistic visualizations of fluid dynamics.

### Dependencies
- C++ compiler
- Blender for visualization
