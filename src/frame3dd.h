/*	FRAME3DD: Static and dynamic structural analysis of 2D & 3D frames and trusses
	Copyright (C) 1992-2008  Henri P. Gavin

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*//**
	@file
	Main functions of the FRAME3DD solver API
*/

#ifndef FRAME_FRAME_H
#define FRAME_FRAME_H

/* note, common.h includes a sneaky #define that rewrites 'float' to 'double'. */
#include "common.h"

/** form the global stiffness matrix */
void assemble_K(
	float **K
	, int DoF, int nM
	, float *x, float *y, float *z, float *r, float *L, float *Le
	, int *J1, int *J2
	, float *Ax, float *Asy, float *Asz
	, float *J, float *Iy, float *Iz, float *E, float *G, float *p
	, int shear, int geom, float **Q
);

/** apply boundary conditions */
void apply_reactions(
	int DoF, int *R
	, float *Dp, float *Fo, float *F, float **K
);

/** solve {F} =   [K]{D} via L D L' decomposition */
void solve_system(
	float **K, float *D, float *F, int DoF, int *ok
);

/** evaluate the member end forces for every member */
void end_forces(
	float **Q, int nM, float *x, float *y, float *z
	, float *L, float *Le
	, int *J1, int *J2, float *Ax, float *Asy, float *Asz
	, float *J, float *Iy, float *Iz, float *E, float *G, float *p
	, float *D, int shear, int geom
);

/** perform an equilibrium check, F returned as reactions */
void equilibrium(	
	float *x, float *y, float *z
	, float *L, int *J1, int *J2, float *F, int *R, float *p
	, float **Q, float **feF, int nM, int DoF, float *err
);

/** assemble global mass matrix from element mass & inertia */
void assemble_M(
	float **M, int DoF, int nJ, int nM,
	float *x, float *y, float *z, float *r, float *L
	, int *J1, int *J2
	, float *Ax, float *J, float *Iy, float *Iz, float *p
	, float *d, float *BMs, float *JMs, float *JMx, float *JMy, float *JMz
	, int lump
);

/** static condensation of stiffness matrix from NxN to nxn */
void condense(
	float **A, int N, int *q, int n, float **Ac
);


/**
	generalized Guyan reduction of mass and stiffness matrices
	matches the response at a particular frequency, sqrt(L)/2/pi
	Guyan, Robert J., "Reduction of Stiffness and Mass Matrices",
	AIAA Journal, Vol. 3, No. 2 (1965) p 380.
*/
void guyan(
	float **M, float **K, int N
	, int *q, int n
	, float **Mc, float **Kc, float w2 
);

/**
	dynamic condensation of mass and stiffness matrices
	matches the response at a set of frequencies

	@NOTE Kc and Mc may be ill-conditioned, and possibly non-positive def.
*/
void dyn_conden(
	float **M, float **K, int N, int *R, int *p, int n
	, float **Mc, float **Kc, float **V, float *f, int *m
);

/**
	release allocated memory
*/
void deallocate( 
	float *x, float *y, float *z, float *r, float *L, float *Le
	, int *J1, int *J2
	, float *Ax, float *Asy, float *Asz, float *J, float *Iy, float *Iz
	, float *E, float *G, float **K, float **Q, float *F, float *D
	, int *R
	, float **W,  float **P,  float **T, float **feF
	, float *Fo,  float *d,  float *BMs, float *JMs
	, float *JMx, float *JMy, float *JMz,  float **M, float *f, float **V
	, int nJ, int nM, int DoF, int modes 
);


#endif /* FRAME_FRAME_H */

