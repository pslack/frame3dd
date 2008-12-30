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
	Frame input/output functions (FRAME3DD native format)
*/

#include "common.h"
#include "microstran/vec3.h"

#include <stdio.h>

/**
	get line into a character string. from K&R.
	@NOTE this is different from the GNU 'getline'.
*/
void frm_getline( FILE *fp, char *s, int lim );


/* re-write input file without comments */
void parse_input(FILE *fp);


/* read input data file			*/
void read_input_data(
	FILE *fp,
	int nJ, int nM, vec3 *pos,
	double *r, double *L, double *Le,
	int *J1, int *J2, 
	double *Ax, double *Asy, double *Asz,
	double *J, double *Iy, double *Iz, double *E, double *G, double *p
);

/** 
	Read data describing run analysis
*/
void read_run_data (
	FILE *fp, int *shear, int *geom,
	char mesh_file[], char plot_file[],
	double *exagg, int *anlyz
);


/**
	Read fixed joint displacement boundary conditions
*/
void read_reaction_data(
	FILE *fp, int DoF, int *nR, int nJ, int *R, int *sumR 
);


/**
	read load information data, form un-restrained load vector
*/
void assemble_loads(
	FILE *fp, 
	int nL, int nJ, vec3 *pos,
	double *L, double *Le, double *Ax, double *Asy, double *Asz, 
	double *Iy, double *Iz, double *E, double *G, 
	double *p, int shear,
	int *J1, int *J2,
	int DoF, 
	int nM, int *nF, int *nW, int *nP, int *nT, int *nD,
	double **Q, double **F_mech, double **F_temp, 
	double ***W, double ***P, double ***T, 
	double ***feF_mech, double ***feF_temp, double **Dp, int *R
);


/**
	read member densities and extra inertial mass data
*/
void read_mass_data(
	FILE *fp, 
	int nJ, int nM, int *nI, 
	double *d, double *BMs, 
	double *JMs, double *JMx, double *JMy, double *JMz, 
	double *L, double *Ax, 
	double *total_mass, double *struct_mass, 
	int *modes, int *Mmethod, int *lump, 
	char modefile[], 
	double *tol, double *shift, int anim[], int *pan 
);


/**
	read matrix condensation information
*/
void read_condense(
	FILE *fp, 
	int nJ, int modes, 
	int *nC, int *Cdof, int *Cmethod, int *q, int *m
);


/**
	write input data to a file
*/
void write_input_data(
	FILE *fp, 
	char *title, int nJ, int nM,  int nL, 
	int *nD, int nR, 
	int *nF, int *nW, int *nP, int *nT,
	vec3 *pos, double *r,
	int *J1, int *J2,
	double *Ax, double *Asy, double *Asz, double *J, double *Iy, double *Iz,
	double *E, double *G, double *p, double **F, double **Dp,
	int *R,
	double ***W, double ***P, double ***T,
	int shear, int anlyz, int geom
);


/**
	save joint displacements and member end forces in a text file	9sep08
*/
void write_static_results(
	FILE *fp, int nJ, int nM, int nL, int lc, int DoF, 
	int *J1, int *J2, 
	double *F, double *D, int *R, 
	double **Q, double err, 
	int ok 
);


/**
	save joint displacements and member end forces in an m-file	9sep08
*/
void write_static_mfile ( char *argv[],
	int nJ, int nM, int nL, int lc, int DoF,
	int *J1, int *J2,
	double *F, double *D, int *R,
	double **Q, double err,
	int ok
);


/**
	save modal frequencies and mode shapes			16aug01
*/
void write_modal_results(
	FILE *fp, 
	int nJ, int nM, int nI, int DoF, 
	double **M, double *f, double **V, 
	double total_mass, double struct_mass, 
	int iter, int sumR, int modes, 
	double shift, int lump, double tol, int ok 
);


/**	
	create mesh data of deformed and undeformed mesh, use gnuplot	22feb99
	useful gnuplot options: set noxtics noytics noztics noborder view nokey
*/
void static_mesh(
	char IO_file[], char meshfile[], char plotfile[], 
	char *title, int nJ, int nM, int nL, int lc, int DoF, 
	vec3 *pos, double *L, 
	int *J1, int *J, double *p, double *D, 
	double exg, int anlyz
);


/**
	create mesh data of the mode-shape meshes, use gnuplot	19oct98
	useful gnuplot options: set noxtics noytics noztics noborder view nokey
*/
void modal_mesh(
	char IO_file[], char meshfile[], char modefile[],
	char plotfile[], char *title,
	int nJ, int nM, int DoF, int modes,
	vec3 *pos, double *L,
	int *J1, int *J2, double *p,
	double **M, double *f, double **V,
	double exg, int anlyz
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
	int nJ, int nM, int DoF, int modes,
	vec3 *pos, double *L, double *p,
	int *J1, int *J2, double *f, double **V,
	double exg,
	int pan
);


/**
	computes cubic deflection functions from beam end deflections
	and beam end rotations.  Saves deflected shapes to a file. 
	These bent shapes are exact for mode-shapes, and for frames
	loaded at their joints.
*/
void bent_beam(
	FILE *fp, int j1, int j2,
	vec3 *pos,
	double L, double p, double *D,
	double exg
);


/** print a set of dots (periods) */
void dots(int n);

