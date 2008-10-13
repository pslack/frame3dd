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

#include "microstran/vec3.h"

/** form the global stiffness matrix */
void assemble_K(
	float **K
	, int DoF, int nM
	, vec3 *pos, float *r, float *L, float *Le
	, int *J1, int *J2
	, float *Ax, float *Asy, float *Asz
	, float *J, float *Iy, float *Iz, float *E, float *G, float *p
	, int shear, int geom, float **Q
);

/** apply boundary conditions */
void apply_reactions(
	int DoF, int *R, float *Dp, float *Fo, float *F, float **K
);

/** solve {F} =   [K]{D} via L D L' decomposition */
void solve_system(
	float **K, float *D, float *F, int DoF, int *ok
);

/** evaluate the member end forces for every member */
void end_forces(
	float **Q, int nM, vec3 *pos
	, float *L, float *Le
	, int *J1, int *J2, float *Ax, float *Asy, float *Asz
	, float *J, float *Iy, float *Iz, float *E, float *G, float *p
	, float *D, int shear, int geom
);

/** perform an equilibrium check, F returned as reactions */
void equilibrium(	
	vec3 *pos
	, float *L, int *J1, int *J2, float *F, int *R, float *p
	, float **Q, float **feF, int nM, int DoF, float *err
);

/** assemble global mass matrix from element mass & inertia */
void assemble_M(
	float **M, int DoF, int nJ, int nM,
	vec3 *pos, float *r, float *L
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
	int nJ, int nM, int nL, int *nF, int *nW, int *nP, int *nT, int DoF,
	int modes,
	vec3 *pos, float *r, float *L, float *Le,
	int *J1, int *J2, int *R,
	float *Ax, float *Asy, float *Asz, float *J, float *Iy, float *Iz,
	float *E, float *G, float *p,
	float ***W, float ***P, float ***T,
	float **Fo_mech, float **Fo_temp, float *Fo_mech_lc, float *Fo_temp_lc,
	float ***feF_mech, float ***feF_temp, float **feF,
	float **Fo, float *Fo_lc, float *F_lc,
	float **K, float **Q,
	float *D, float *dD, float *Dp,
	float *d, float *BMs, float *JMs, float *JMx, float *JMy, float *JMz,
	float **M, float *f, float **V, 
	int *q, int *m
);


#endif /* FRAME_FRAME_H */

