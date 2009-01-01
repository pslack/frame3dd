/*
 FRAME3DD:
 Static and dynamic structural analysis of 2D and 3D frames and trusses with
 elastic and geometric stiffness.
 ---------------------------------------------------------------------------
 http://www.duke.edu/~hpgavin/frame/
 ---------------------------------------------------------------------------
 Copyright (C) 1992-2009  Henri P. Gavin
 
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
	double **K,		/**< stiffness matrix		*/
	int DoF, int nB,	/**< degrees of freedom, number of beams */
	vec3 *xyz,		/**< XYZ locations of every joint	*/
	float *r,		/**< rigid radius of every joint	*/
	double *L, double *Le,	/**< length of each beam element, effective */
	int *J1, int *J2,	/**< joint connectivity			*/
	float *Ax, float *Asy, float *Asz,	/**< section areas	*/
	float *J, float *Iy, float *Iz,	/**< section inertias	*/
	float *E, float *G,	/**< elastic and shear moduli		*/
	float *p,		/**< roll angle, radians		*/
	int shear,		/**< 1: include shear deformation, 0: don't */
	int geom,		/**< 1: include goemetric stiffness, 0: don't */
	double **Q		/**< beam element end forces for every beam */
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
	int *ok		/**< indicates positive definite stiffness matrix */
);


/** evaluate the member end forces for every member */
void end_forces(
	double **Q,	/**< beam element end forces for every beam	*/
	int nB,		/**< number of beam elements			*/
	vec3 *xyz,	/** XYZ locations of each joint			*/
	double *L, double *Le,	/**< length of each beam element, effective */
	int *J1, int *J2,	/**< joint connectivity			*/
	float *Ax, float *Asy, float *Asz,	/**< section areas	*/
	float *J, float *Iy, float *Iz,	/**< section area inertias	*/
	float *E, float *G,	/**< elastic and shear moduli		*/
	float *p,		/**< roll angle, radians		*/
	double *D,	/**< displacement vector			*/
	int shear,	/**< 1: include shear deformation, 0: don't */
	int geom	/**< 1: include goemetric stiffness, 0: don't */
);


/** perform an equilibrium check, F returned as reactions */
void equilibrium(	
	vec3 *xyz,	/** XYZ locations of each joint			*/
	double *L,	/**< length of each beam element, effective	*/
	int *J1, int *J2, /**< joint connectivity			*/
	double *F,	/**< load vector				*/
	int *R,		/**< R[i]=1: DoF i is fixed, R[i]=0: DoF i is free */
	float *p,	/**< roll angle, radians			*/
	double **Q,	/**< beam element end forces for every beam	*/
	double **feF,	/**< fixed end forces for every beam element	*/
	int nB,		/**< number of beam elements			*/
	int DoF,	/**< number of degrees of freedom		*/
	double *err	/**< root mean squared equilibrium error	*/
);


/** assemble global mass matrix from element mass & inertia */
void assemble_M(
	double **M,	/**< mass matrix				*/
	int DoF,	/**< number of degrees of freedom		*/
	int nJ, int nB,	/**< number of joints, number of beam elements	*/
	vec3 *xyz,	/** XYZ locations of each joint			*/
	float *r,	/**< rigid radius of every joint		*/
	double *L,	/**< length of each beam element, effective	*/
	int *J1, int *J2, /**< joint connectivity			*/
	float *Ax,	/**< joint connectivity				*/
	float *J, float *Iy, float *Iz,	/**< section area inertias	*/
	float *p,	/**< roll angle, radians			*/
	float *d,	/**< beam element density			*/
	float *BMs,	/**< extra beam mass				*/
	float *JMs,	/**< joint mass					*/
	float *JMx, float *JMy, float *JMz,	/**< joint inertias	*/
	int lump	/**< 1: lumped mass matrix, 0: consistent mass	*/
);


/** static condensation of stiffness matrix from NxN to nxn */
void condense(
	double **A,	/**< a square matrix				*/
	int N,		/**< the dimension of the matrix		*/
	int *q,		/**< list of matrix indices to retain		*/
	int n,		/**< the dimension of the condensed matrix	*/
	double **Ac	/**< the condensed matrix			*/
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
	double w2 		/**< matched value of frequency squared	*/
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
	int *m		/**< list of modes to match in the condensed model */
);


/**
	release allocated memory
*/
void deallocate( 
	int nJ, int nB, int nL, int *nF, int *nW, int *nP, int *nT, int DoF,
	int modes,
	vec3 *xyz, float *r, double *L, double *Le,
	int *J1, int *J2, int *R,
	float *Ax, float *Asy, float *Asz,
	float *J, float *Iy, float *Iz,
	float *E, float *G,
	float *p,
	float ***W, float ***P, float ***T,
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

