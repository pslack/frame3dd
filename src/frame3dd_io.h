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
	Frame input/output functions (FRAME3DD native format)
*/

#include "common.h"
#include "time.h"
#include "microstran/vec3.h"

#include <stdio.h>

/**
	get line into a character string. from K&R.
	@NOTE this is different from the GNU 'getline'.
*/
void frame3dd_getline( FILE *fp, char *s, int lim );


/* re-write input file without comments or commas */
void parse_input(FILE *fp);


/* read input data file			*/
void read_input_data(
	FILE *fp,	/**< input data file pointer			*/
	int nJ, int nB,	/**< number of joints, number of beam elements	*/
	vec3 *xyz,	/**< XYZ coordinates of each joint		*/
	float *r,	/**< rigid radius of each joint			*/
	double *L, double *Le,	/**< length of each beam element, effective */
	int *J1, int *J2, 	/**< joint connectivity			*/
	float *Ax, float *Asy, float *Asz,	/**< section areas	*/
	float *J, float *Iy, float *Iz,	/**< section inertias	*/
	float *E, float *G,	/**< elastic moduli and shear moduli	*/
	float *p		/**< roll angle of each beam (radians)	*/
);

/** 
	Read data controlling certain aspects of the analysis
*/
void read_run_data (
	FILE *fp,	/**< input data file pointer			*/
	int *shear,	/**< 1: include shear deformations, 0: don't	*/
	int *geom,	/**< 1: include geometric stiffness, 0: don't	*/
	char mesh_file[],	/**< file name for mesh data output	*/
	char plot_file[],	/**< file name for Gnuplot script	*/
	double *exagg,		/**< factor for deformation exaggeration */
	int *anlyz	/* 1: perform elastic analysis, 0: don't	*/
);


/**
	Read fixed joint displacement boundary conditions
*/
void read_reaction_data(
	FILE *fp,	/**< input data file pointer			*/
	int DoF,	/**< number of degrees of freedom		*/
	int nJ,		/**< number of joints				*/
	int *nR,	/**< number of joints with reactions		*/
	int *R,		/**< R[i]=1: DoF i is fixed, R[i]=0: DoF i is free */
	int *sumR 	/**< sum of vector R				*/
);


/**
	read load information data, form un-restrained load vector
*/
void read_and_assemble_loads(
	FILE *fp,	/**< input data file pointer			*/
	int nJ,		/**< number of joints				*/
	int nB,		/**< number of beam elements			*/
	int nL,		/**< number of load cases	i		*/
	int DoF,	/**< number of degrees of freedom		*/
	vec3 *xyz,	/**< XYZ coordinates of each joint		*/
	double *L, double *Le,	/**< length of each beam element, effective */
	int *J1, int *J2, 	/**< joint connectivity			*/
	float *Ax, float *Asy, float *Asz,	/**< section areas	*/
	float *Iy, float *Iz,	/**< section inertias			*/
	float *E, float *G,	/**< elastic moduli and shear moduli	*/
	float *p,		/**< roll angle of each beam (radians)	*/
	int *R,		/**< R[i]=1: DoF i is fixed, R[i]=0: DoF i is free */
	int shear,	/**< 1: include shear deformations, 0: don't	*/
	int *nF, int *nW, int *nP, int *nT, int *nD,	/**< number of loads */
	double **Q,		/**< beam element end forces, every beam */
	double **F_mech, 	/**< mechanical loads			*/
	double **F_temp, 	/**< thermal loads			*/
	float ***W, float ***P, float ***T, 	/**< loads		*/
	float **Dp,		/**< prescribed displacements at rctns	*/
	double ***feF_mech, double ***feF_temp	/**< fixed end forces	*/
);


/**
	read member densities and extra inertial mass data
*/
void read_mass_data(
	FILE *fp,	/**< input data file pointer			*/
	int nJ, int nB,	/**< number of joints, number of beams		*/
	int *nI,	/**< number of beams with extra inertia		*/
	float *d, float *BMs, /**< beam density, extra beam mass	*/
	float *JMs, float *JMx, float *JMy, float *JMz, /**< joint inertia*/
	double *L,	/**< length of each beam element		*/
	float *Ax, 	/**< cross section area of each beam element	*/
	double *total_mass,	/**< total mass of structure and extra mass */
	double *struct_mass, 	/**< mass of structural elements	*/
	int *nM,	/**< number of modes to find			*/
	int *Mmethod, 	/**< modal analysis method			*/
	int *lump,	/**< 1: use lumped mass matrix, 0: consistent mass */
	char modefile[], /**< filename for mode shape data for plotting	*/	
	double *tol,	/**< convergence tolerance for mode shapes	*/
	double *shift,	/**< frequency shift for unrestrained frames	*/
	int *anim,	/**< list of modes to be graphically animated	*/
	int *pan 	/**< 1: pan viewpoint during animation, 0: don't */
);


/**
	read matrix condensation information
*/
void read_condensation_data(
	FILE *fp,	/**< input data file pointer			*/
	int nJ, int nM, 	/**< number of joints, number of modes	*/
	int *nC,	/**< number of joints with condensed DoF's	*/
	int *Cdof,	/**< list of DoF's retained in condensed model	*/
	int *Cmethod,	/**< matrix conden'n method, static, Guyan, dynamic*/
	int *q,		/**< list of retained degrees of freedom	*/
	int *m		/**< list of retained modes in dynamic condensation */
);


/**
	write input data to a file
*/
void write_input_data(
	FILE *fp,	/**< input data file pointer			*/
	char *title, int nJ, int nB,  int nL, 
	int *nD, int nR, 
	int *nF, int *nW, int *nP, int *nT,
	vec3 *xyz, float *r,
	int *J1, int *J2,
	float *Ax, float *Asy, float *Asz,
	float *J, float *Iy, float *Iz,
	float *E, float *G, float *p,
	double **F, float **Dp,
	int *R,
	float ***W, float ***P, float ***T,
	int shear, int anlyz, int geom
);


/**
	save joint displacements and member end forces in a text file	9sep08
*/
void write_static_results(
	FILE *fp,
	int nJ, int nB, int nL, int lc, int DoF, 
	int *J1, int *J2, 
	double *F, double *D, int *R, double **Q,
	double err, int ok 
);


/**
	save joint displacements and member end forces in a .CSV file   31dec08
*/
void write_static_csv(
	char *argv[], char *title,
	int nJ, int nB, int nL, int lc, int DoF, 
	int *J1, int *J2, 
	double *F, double *D, int *R, double **Q,
	double err, int ok 
);


/**
	save joint displacements and member end forces in an m-file	9sep08
*/
void write_static_mfile ( 
	char *argv[], char *title,
	int nJ, int nB, int nL, int lc, int DoF,
	int *J1, int *J2,
	double *F, double *D, int *R, double **Q,
	double err, int ok
);


/**
	save modal frequencies and mode shapes			16aug01
*/
void write_modal_results(
	FILE *fp, 
	int nJ, int nB, int nI, int DoF, 
	double **M, double *f, double **V, 
	double total_mass, double struct_mass, 
	int iter, int sumR, int nM, 
	double shift, int lump, double tol, int ok 
);


/**	
	create mesh data of deformed and undeformed mesh, use gnuplot	22feb99
	useful gnuplot options: set noxtics noytics noztics noborder view nokey
*/
void static_mesh(
	char IO_file[], char meshfile[], char plotfile[], 
	char *title, int nJ, int nB, int nL, int lc, int DoF, 
	vec3 *xyz, double *L, 
	int *J1, int *J, float *p, double *D, 
	double exagg, int anlyz
);


/**
	create mesh data of the mode-shape meshes, use gnuplot	19oct98
	useful gnuplot options: set noxtics noytics noztics noborder view nokey
*/
void modal_mesh(
	char IO_file[], char meshfile[], char modefile[],
	char plotfile[], char *title,
	int nJ, int nB, int DoF, int nM,
	vec3 *xyz, double *L,
	int *J1, int *J2, float *p,
	double **M, double *f, double **V,
	double exagg, int anlyz
);


/**
	create mesh data of animated mode-shape meshes, use gnuplot	16dec98
	useful gnuplot options: set noxtics noytics noztics noborder view nokey
	mpeg movie example:   % convert mesh_file-03-f-*.ps mode-03.mpeg
	... requires ImageMagick and mpeg2vidcodec packages
*/
void animate(
	char IO_file[], char meshfile[], char modefile[], char plotfile[],
	char *title,
	int anim[],
	int nJ, int nB, int DoF, int nM,
	vec3 *xyz, double *L, float *p,
	int *J1, int *J2, double *f, double **V,
	double exagg, int pan
);


/**
	computes cubic deflection functions from beam end deflections
	and beam end rotations.  Saves deflected shapes to a file. 
	These bent shapes are exact for mode-shapes, and for frames
	loaded at their joints.
*/
void bent_beam(
	FILE *fp, int j1, int j2,
	vec3 *xyz,
	double L, float p, double *D,
	double exagg
);


/**
	return the file extension, including the period (.)
                return 1 if the extension is ".csv"
                return 2 if the extension is ".fmm"
                return 0 otherwise

*/
int get_file_ext( char *filename, char *ext );



/** print a set of dots (periods) */
void dots( int n );

