/*
 FRAME3DD:
 Static and dynamic structural analysis of 2D and 3D frames and trusses with
 elastic and geometric stiffness.
 ---------------------------------------------------------------------------
 http://www.duke.edu/~hpgavin/frame/
 ---------------------------------------------------------------------------
 Copyright (C) 1992-2008  Henri P. Gavin
 
    FRAME3DD is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FRAME3DD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FRAME3DD.  If not, see <http://www.gnu.org/licenses/>.
*//**
	@file
	Main functions of the FRAME3DD solver API
*/

#ifndef FRAME_FRAME_H
#define FRAME_FRAME_H

/* for Micro-Stran compatability, structure for cartesian vectors */
#include "microstran/vec3.h"

/* maximum number of load cases */
#define _NL_ 128


/** form the global stiffness matrix */
void assemble_K(
	double **K,
	int DoF, int nM,
	vec3 *pos, double *r, double *L, double *Le,
	int *J1, int *J2,
	double *Ax, double *Asy, double *Asz,
	double *J, double *Iy, double *Iz, double *E, double *G, double *p,
	int shear, int geom, double **Q
);


/** apply boundary conditions */
void apply_reactions(
	int DoF, int *R, double **Dp, int lc, 
	double *Fo, double *F, double **K, char tm
);


/** solve {F} =   [K]{D} via L D L' decomposition */
void solve_system(
	double **K, double *D, double *F, int DoF, int *ok
);


/** evaluate the member end forces for every member */
void end_forces(
	double **Q, int nM, vec3 *pos,
	double *L, double *Le,
	int *J1, int *J2, double *Ax, double *Asy, double *Asz,
	double *J, double *Iy, double *Iz, double *E, double *G, double *p,
	double *D, int shear, int geom
);


/** perform an equilibrium check, F returned as reactions */
void equilibrium(	
	vec3 *pos,
	double *L, int *J1, int *J2, double *F, int *R, double *p,
	double **Q, double **feF, int nM, int DoF, double *err
);


/** assemble global mass matrix from element mass & inertia */
void assemble_M(
	double **M, int DoF, int nJ, int nM,
	vec3 *pos, double *r, double *L,
	int *J1, int *J2,
	double *Ax, double *J, double *Iy, double *Iz, double *p,
	double *d, double *BMs, double *JMs, double *JMx, double *JMy, double *JMz,
	int lump
);


/** static condensation of stiffness matrix from NxN to nxn */
void condense(
	double **A, int N, int *q, int n, double **Ac
);


/**
	generalized Guyan reduction of mass and stiffness matrices
	matches the response at a particular frequency, sqrt(L)/2/pi
	Guyan, Robert J., "Reduction of Stiffness and Mass Matrices",
	AIAA Journal, Vol. 3, No. 2 (1965) p 380.
*/
void guyan(
	double **M, double **K, int N,
	int *q, int n,
	double **Mc, double **Kc, double w2 
);


/**
	dynamic condensation of mass and stiffness matrices
	matches the response at a set of frequencies

	@NOTE Kc and Mc may be ill-conditioned, and possibly non-positive def.
*/
void dyn_conden(
	double **M, double **K, int N, int *R, int *p, int n,
	double **Mc, double **Kc, double **V, double *f, int *m
);


/**
	release allocated memory
*/
void deallocate( 
	int nJ, int nM, int nL, int *nF, int *nW, int *nP, int *nT, int DoF,
	int modes,
	vec3 *pos, double *r, double *L, double *Le,
	int *J1, int *J2, int *R,
	double *Ax, double *Asy, double *Asz, double *J, double *Iy, double *Iz,
	double *E, double *G, double *p,
	double ***W, double ***P, double ***T,
	double **Fo_mech, double **Fo_temp, double *Fo_mech_lc, double *Fo_temp_lc,
	double ***feF_mech, double ***feF_temp, double **feF,
	double **Fo, double *Fo_lc, double *F_lc,
	double **K, double **Q,
	double *D, double *dD, double **Dp,
	double *d, double *BMs, double *JMs, double *JMx, double *JMy, double *JMz,
	double **M, double *f, double **V, 
	int *q, int *m
);


#endif /* FRAME_FRAME_H */

