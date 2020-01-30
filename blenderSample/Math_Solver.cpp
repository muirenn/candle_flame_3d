#include "stdafx.h"
#include "Math_Solver.h"

/*
Velocity equation: 
∂u/∂t=-(u∙∇)u-1/ρ ∇p+ν∇^2 u+f
*/



void lin_solve(property prop, double* x, double* x0, double a, double c, int lin_itr, unsigned long int n) {
	unsigned long int i, j, k, l;
	for (l = 0; l < lin_itr; l++)
	{
		for (i = 1; i <= n; i++)
		{
			for (j = 1; j <= n; j++)
			{
				for (k = 1; k <= n; k++)
				{
					x[IDX(i, j, k, n)] = (x0[IDX(i, j, k, n)] +
						a * (
							x[IDX(i - 1, j, k, n)] +
							x[IDX(i + 1, j, k, n)] +
							x[IDX(i, j - 1, k, n)] +
							x[IDX(i, j + 1, k, n)] +
							x[IDX(i, j, k - 1, n)] +
							x[IDX(i, j, k + 1, n)]
							)
						) / c;
				}
			}
		}
		set_bnd(prop, x, n);
	}
}

void set_bnd(property prop, double* x, unsigned long int n)
{
	unsigned long int i, j, k;
	for (i = 1; i <= n; i++)
	{
		for (j = 1; j <= n; j++)
		{
			x[IDX(i, j, 0, n)] = (prop == 3 ? -x[IDX(i, j, 1, n)] : x[IDX(i, j, 1, n)]);
			x[IDX(i, j, n + 1, n)] = (prop == 3 ? -x[IDX(i, j, n, n)] : x[IDX(i, j, n, n)]);
		}
	}
	for (i = 1; i <= n; i++)
	{
		for (k = 1; k <= n; k++)
		{
			x[IDX(i, 0, k, n)] = (prop == 2 ? -x[IDX(i, 1, k, n)] : x[IDX(i, 1, k, n)]);
			x[IDX(i, n + 1, k, n)] = (prop == 2 ? -x[IDX(i, n, k, n)] : x[IDX(i, n, k, n)]);
		}
	}
	for (j = 1; j <= n; j++)
	{
		for (k = 1; k <= n; k++)
		{
			x[IDX(0, j, k, n)] = (prop == 1 ? -x[IDX(1, j, k, n)] : x[IDX(1, j, k, n)]);
			x[IDX(n + 1, j, k, n)] = (prop == 1 ? -x[IDX(n, j, k, n)] : x[IDX(n, j, k, n)]);
		}
	}

	x[IDX(0, 0, 0, n)] = 0.33 * (
		x[IDX(1, 0, 0, n)] +
		x[IDX(0, 1, 0, n)] +
		x[IDX(0, 0, 1, n)]);
	x[IDX(0, n + 1, 0, n)] = 0.33 * (
		x[IDX(1, n + 1, 0, n)] +
		x[IDX(0, n, 0, n)] +
		x[IDX(0, n + 1, 1, n)]);
	x[IDX(0, 0, n + 1, n)] = 0.33 * (
		x[IDX(1, 0, n + 1, n)] +
		x[IDX(0, 1, n + 1, n)] +
		x[IDX(0, 0, n, n)]);
	x[IDX(0, n + 1, n + 1, n)] = 0.33 * (
		x[IDX(1, n + 1, n + 1, n)] +
		x[IDX(0, n, n + 1, n)] +
		x[IDX(0, n + 1, n, n)]);
	x[IDX(n + 1, 0, 0, n)] = 0.33 * (
		x[IDX(n, 0, 0, n)] +
		x[IDX(n + 1, 1, 0, n)] +
		x[IDX(n + 1, 0, 1, n)]);
	x[IDX(n + 1, n + 1, 0, n)] = 0.33 * (
		x[IDX(n, n + 1, 0, n)] +
		x[IDX(n + 1, n, 0, n)] +
		x[IDX(n + 1, n + 1, 1, n)]);
	x[IDX(n + 1, 0, n + 1, n)] = 0.33 * (
		x[IDX(n, 0, n + 1, n)] +
		x[IDX(n + 1, 1, n + 1, n)] +
		x[IDX(n + 1, 0, n, n)]);
	x[IDX(n + 1, n + 1, n + 1, n)] = 0.33 * (
		x[IDX(n, n + 1, n + 1, n)] +
		x[IDX(n + 1, n, n + 1, n)] +
		x[IDX(n + 1, n + 1, n, n)]);
}

// For mass-concerving property of velocity; Poisson equation solver
void project(double* u, double* v, double* w, double* p, double* div, int lin_itr, unsigned long int n)
{
	unsigned long int i, j, k;
	//double* p = new double[52*52*52];
	//double* div = new double[52 * 52 * 52];
	double h = 1.0 / n;

	for (i = 1; i <= n; i++)
	{
		for (j = 1; j <= n; j++)
		{
			for (k = 1; k <= n; k++)
			{
				div[IDX(i, j, k, n)] =
					-0.5f * h * (
						u[IDX(i + 1, j, k, n)] -
						u[IDX(i - 1, j, k, n)] +
						v[IDX(i, j + 1, k, n)] -
						v[IDX(i, j - 1, k, n)] +
						w[IDX(i, j, k + 1, n)] -
						w[IDX(i, j, k - 1, n)]);
				p[IDX(i, j, k, n)] = 0;
			}
		}
	}
	// p - gradient field?
	set_bnd(density, div, n);
	set_bnd(density, p, n);
	lin_solve(density, p, div, 1, 6, lin_itr, n);

	for (i = 1; i <= n; i++)
	{
		for (j = 1; j <= n; j++)
		{
			for (k = 1; k <= n; k++)
			{
				u[IDX(i, j, k, n)] -= 0.5f * (
					p[IDX(i + 1, j, k, n)] -
					p[IDX(i - 1, j, k, n)]
					) * n;
				v[IDX(i, j, k, n)] -= 0.5f * (
					p[IDX(i, j + 1, k, n)] -
					p[IDX(i, j - 1, k, n)]
					) * n;
				w[IDX(i, j, k, n)] -= 0.5f * (
					p[IDX(i, j, k + 1, n)] -
					p[IDX(i, j, k - 1, n)]
					) * n;
			}
		}
	}
	set_bnd(velocity_x, u, n);
	set_bnd(velocity_y, v, n);
	set_bnd(velocity_z, w, n);

	//delete[] div;
	//delete[] p;
}


// + ν∇^2 x
void diffuse(double dt, double diff, property prop, double* x, double* x0, int linitr, unsigned long int n)
{
	double a = dt * diff * n * n * n;
	lin_solve(prop, x, x0, a, 1 + 6 * a, linitr, n);
}


// - (u∙∇)d
void advect(double dt, property prop, double* d, double* d0, double* u, double* v, double* w, int n)
{
	unsigned long int i, j, k;
	double x, y, z, s0, t0, r0, s1, t1, r1, dt0, i0, j0, k0, i1, j1, k1;

	dt0 = dt * n;

	for (i = 1; i <= n; i++)
	{
		for (j = 1; j <= n; j++)
		{
			for (k = 1; k <= n; k++)
			{
				x = static_cast<double>(i) - dt0 * u[IDX(i, j, k, n)];
				y = static_cast<double>(j) - dt0 * v[IDX(i, j, k, n)];
				z = static_cast<double>(k) - dt0 * w[IDX(i, j, k, n)];
				if (x < 0.5) x = 0.5;
				if (x > n + 0.5) x = static_cast<double>(n) + 0.5;
				i0 = floor(x);
				i1 = i0 + 1;
				if (y < 0.5) y = 0.5;
				if (y > n + 0.5) y = static_cast<double>(n) + 0.5;
				j0 = floor(y);
				j1 = j0 + 1;
				if (z < 0.5) z = 0.5;
				if (z > n + 0.5) z = static_cast<double>(n) + 0.5;
				k0 = floor(z);
				k1 = k0 + 1;

				s1 = x - i0;
				s0 = 1 - s1;
				t1 = y - j0;
				t0 = 1 - t1;
				r1 = z - k0;
				r0 = 1 - r1;

				int i0i = static_cast<int>(i0);
				int i1i = static_cast<int>(i1);
				int j0i = static_cast<int>(j0);
				int j1i = static_cast<int>(j1);
				int k0i = static_cast<int>(k0);
				int k1i = static_cast<int>(k1);

				const double res =
					s0 * t0 * r0 * d0[IDX(i0i, j0i, k0i, n)] +
					s0 * t0 * r1 * d0[IDX(i0i, j0i, k1i, n)] +
					s0 * t1 * r0 * d0[IDX(i0i, j1i, k0i, n)] +
					s0 * t1 * r1 * d0[IDX(i0i, j1i, k1i, n)] +
					s1 * t0 * r0 * d0[IDX(i1i, j0i, k0i, n)] +
					s1 * t0 * r1 * d0[IDX(i1i, j0i, k1i, n)] +
					s1 * t1 * r0 * d0[IDX(i1i, j1i, k0i, n)] +
					s1 * t1 * r1 * d0[IDX(i1i, j1i, k1i, n)];

				d[IDX(i, j, k, n)] = res;
			}
		}
	}

	set_bnd(prop, d, n);
}