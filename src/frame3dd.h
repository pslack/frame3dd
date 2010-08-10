/*
 This file is part of FRAME3DD:
 Static and dynamic structural analysis of 2D and 3D frames and trusses with
 elastic and geometric stiffness.
 ---------------------------------------------------------------------------
 http://frame3dd.sourceforge.net/
 ---------------------------------------------------------------------------
 Copyright (C) 1992-2010  Henri P. Gavin
 
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
#define _NL_ 32


/** form the global stiffness matrix */
void assemble_K(
	double **K,		/**< stiffness matrix			*/
	int DoF,		/**< number of degrees of freedom	*/
	int nE,			/**< number of frame elements		*/
	vec3 *xyz,		/**< XYZ locations of every joint	*/
	float *r,		/**< rigid radius of every joint	*/
	double *L, double *Le,	/**< length of each frame element, effective */
	int *J1, int *J2,	/**< joint connectivity			*/
	float *Ax, float *Asy, float *Asz,	/**< section areas	*/
	float *Jx, float *Iy, float *Iz,	/**< section inertias	*/
	float *E, float *G,	/**< elastic and shear moduli		*/
	float *p,		/**< roll angle, radians		*/
	int shear,		/**< 1: include shear deformation, 0: don't */
	int geom,		/**< 1: include goemetric stiffness, 0: don't */
	double **Q,		/**< frame element end forces		*/
	int debug		/**< 1: write element stiffness matrices*/
);


/** apply boundary conditions */
void apply_reactions(
	int DoF,	/**< number of degrees of freedom		*/
	int *R,		/**< R[i]=1: DoF i is fixed, R[i]=0: DoF i is free */
	float *Dp,	/**< prescribed displacements, each DoF		*/	
	double *Fo,	/**< fixed end forces for unrestrained frame	*/
	double *F,	/**< load vector for restrained frame		*/
	double **K,	/**< stiffness matrix				*/
	char tm		/**< tm='t': thermal loads; tm='m': mech. loads	*/
);


/** solve {F} =   [K]{D} via L D L' decomposition */
void solve_system(
	double **K,	/**< stiffness matrix for the restrained frame	*/
	double *D,	/**< displacement vector to be solved		*/
	double *F,	/**< load vector				*/
	int DoF,	/**< number of degrees of freedom		*/
	int *ok,	/**< indicates positive definite stiffness matrix */
	int verbose	/**< 1: copious screen output; 0: none		*/
);


/** evaluate the member end forces for every member */
void end_forces(
	double **Q,	/**< frame element end forces			*/
	int nE,		/**< number of frame elements			*/
	vec3 *xyz,	/** XYZ locations of each joint			*/
	double *L, double *Le,	/**< length of each frame element, effective */
	int *J1, int *J2,	/**< joint connectivity			*/
	float *Ax, float *Asy, float *Asz,	/**< section areas	*/
	float *Jx, float *Iy, float *Iz,	/**< section area inertias */
	float *E, float *G,	/**< elastic and shear moduli		*/
	float *p,		/**< roll angle, radians		*/
	double *D,	/**< displacement vector			*/
	int shear,	/**< 1: include shear deformation, 0: don't */
	int geom	/**< 1: include goemetric stiffness, 0: don't */
);


/** perform an equilibrium check, F returned as reactions */
void equilibrium(	
	vec3 *xyz,	/**< XYZ locations of each joint		*/
	double *L,	/**< length of each frame element, effective	*/
	int *J1, int *J2, /**< joint connectivity			*/
	double *F,	/**< load vector				*/
	int *R,		/**< R[i]=1: DoF i is fixed, R[i]=0: DoF i is free */
	float *p,	/**< roll angle, radians			*/
	double **Q,	/**< frame element end forces			*/
	double **feF,	/**< fixed end forces for every frame element	*/
	int nE,		/**< number of frame elements			*/
	int DoF,	/**< number of degrees of freedom		*/
	double *err,	/**< root mean squared equilibrium error	*/
	int verbose	/**< 1: copious screen output; 0: none		*/
);


/** assemble global mass matrix from element mass & inertia */
void assemble_M(
	double **M,	/**< mass matrix				*/
	int DoF,	/**< number of degrees of freedom		*/
	int nJ, int nE,	/**< number of joints, number of frame elements	*/
	vec3 *xyz,	/** XYZ locations of each joint			*/
	float *r,	/**< rigid radius of every joint		*/
	double *L,	/**< length of each frame element, effective	*/
	int *J1, int *J2, /**< joint connectivity			*/
	float *Ax,	/**< joint connectivity				*/
	float *Jx, float *Iy, float *Iz,	/**< section area inertias*/
	float *p,	/**< roll angle, radians			*/
	float *d,	/**< frame element density			*/
	float *BMs,	/**< extra frame element mass			*/
	float *JMs,	/**< joint mass					*/
	float *JMx, float *JMy, float *JMz,	/**< joint inertias	*/
	int lump,	/**< 1: lumped mass matrix, 0: consistent mass	*/
	int debug	/**< 1: write element mass matrices	 	*/
);


/** static condensation of stiffness matrix from NxN to nxn */
void condense(
	double **A,	/**< a square matrix				*/
	int N,		/**< the dimension of the matrix		*/
	int *q,		/**< list of matrix indices to retain		*/
	int n,		/**< the dimension of the condensed matrix	*/
	double **Ac,	/**< the condensed matrix			*/
	int verbose	/**< 1: copious screen output; 0: none		*/
);


/**
	generalized Guyan reduction of mass and stiffness matrices
	matches the response at a particular frequency, sqrt(L)/2/pi
	Guyan, Robert J., "Reduction of Stiffness and Mass Matrices",
	AIAA Journal, Vol. 3, No. 2 (1965) p 380.
*/
void guyan(
	double **M, double **K,	/**< mass and stiffness matrices	*/
	int N,			/**< dimension of the matrices, DoF	*/
	int *q,			/**< list of degrees of freedom to retain */
	int n,			/**< dimension of the condensed matrices */
	double **Mc, double **Kc,	/**< the condensed matrices	*/
	double w2,		/**< matched value of frequency squared	*/
	int verbose	/**< 1: copious screen output; 0: none		*/
);


/**
	dynamic condensation of mass and stiffness matrices
	matches the response at a set of frequencies

	@NOTE Kc and Mc may be ill-conditioned, and xyzsibly non-positive def.
*/
void dyn_conden(
	double **M, double **K,	/**< mass and stiffness matrices	*/
	int N,			/**< dimension of the matrices, DoF	*/
	int *R,		/**< R[i]=1: DoF i is fixed, R[i]=0: DoF i is free */
	int *p,		/**< list of primary degrees of freedom		*/
	int n,		/**< the dimension of the condensed matrix	*/
	double **Mc, double **Kc,	/**< the condensed matrices	*/
	double **V, double *f,	/**< mode shapes and natural frequencies*/
	int *m,		/**< list of modes to match in the condensed model */
	int verbose	/**< 1: copious screen output; 0: none		*/
);


/**
	release allocated memory
*/
void deallocate( 
	int nJ, int nE, int nL, int *nF, int *nU, int *nW, int *nP, int *nT, int DoF,
	int modes,
	vec3 *xyz, float *r, double *L, double *Le,
	int *J1, int *J2, int *R,
	float *Ax, float *Asy, float *Asz,
	float *Jx, float *Iy, float *Iz,
	float *E, float *G,
	float *p,
	float ***U, float ***W, float ***P, float ***T,
	float **Dp,
	double **Fo_mech, double **Fo_temp,
	double ***feF_mech, double ***feF_temp, double **feF,
	double **Fo, double *F_lc,
	double **K, double **Q,
	double *D, double *dD,
	float *d, float *BMs,
	float *JMs, float *JMx, float *JMy, float *JMz,
	double **M, double *f, double **V, 
	int *q, int *m
);

/**
	relative norm
	compute the relative 2-norm between two vectors N and D
	return  ( sqrt(sum(N[i]*N[i]) / sqrt(D[i]*D[i]) )
 */
double rel_norm( double *N, double *D, int n );


#endif /* FRAME_FRAME_H */

