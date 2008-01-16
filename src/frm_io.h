/*	FRAME: Static and dynamic structural analysis of 2D & 3D frames and trusses
	Copyright (C) 1992-2007  Henri P. Gavin

	This program is free software; you may redistribute it and/or
	modify it under the terms of the GNU General Public License
	as published by the Free Software Foundation; either version 2
	of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "common.h"

#include <stdio.h>

/**
	get line into a character string. from K&R.
	@NOTE this is different from the GNU 'getline'.
*/
void frm_getline(
		FILE	*fp
		, char    *s
		, int     lim
);

/* read input data file			*/
void read_input(
		FILE *fp
		, int nJ, int nM, float *x, float *y, float *z
		, float *r, float *L, float *Le
		, int *J1, int *J2, int *anlyz, int *geom, float **Q
		, float *Ax, float *Asy, float *Asz
		, float *J, float *Iy, float *Iz, float *E, float *G, float *p
		, int *shear, char meshfile[], char plotfile[], float *exagg
);

/* re-write input file without comments */
void parse_input(FILE *fp);

/**
	read load information data, form un-restrained load vector
*/
void read_loads(
		FILE *fp
		, int nJ, float *x, float *y, float *z
		, float *L, float *Le, float *Ax, float *Asy, float *Asz
		, float *Iy, float *Iz, float *E, float *G
		, float *p, int shear
		, int *J1, int *J2
		, int DoF
		, int nM, int *nF, int *nW, int *nP, int *nT
		, float *F_mech, float *F_temp
		, float **W, float **P, float **T, float **feF_mech, float **feF_temp
);

/**
	Read fixed joint displacement boundary conditions
*/
void read_reactions(
		FILE *fp
		, int DoF, int *nD, int *nR
		, int nJ, float *Dp, int *R, int *sumR 
);

/**
	read member densities and extra inertial mass data
*/
void read_masses(
		FILE *fp
		, int nJ, int nM, int *nI
		, float *d, float *BMs, float *JMs, float *JMx, float *JMy, float *JMz
		, float *L, float *Ax
		, float *total_mass, float *struct_mass
		, int *modes, int *Mmethod, int *lump
		, char modefile[]
		, float *tol, float *shift, int anim[], int *pan
);

/**
	read matrix condensation information
*/
void read_condense(
		FILE *fp
		, int nJ, int modes
		, int *nC, int *Cdof, int *Cmethod, int *q, int *m
);

/**
	save input data
*/
void control_data(
		FILE *fp
		, char *title, int nJ, int nM, int nF, int nD, int nR, int nW, int nP,int nT
		, float *x, float *y, float *z, float *r
		, int *J1, int *J2
		, float *Ax, float *Asy, float *Asz, float *J, float *Iy, float *Iz
		, float *E, float *G, float *p, float *F, float *Dp
		, int *R
		, float **W, float **P, float **T
		, int shear, int anlyz, int geom
);

/**
	save joint displacements and member end forces		15oct98
*/
void save_results(
		FILE *fp, int nJ, int nM, int DoF, int *J1, int *J2
		, float *F, float *D, int *R
		, float **Q, float err
		, int ok 
);

/**
	save modal frequencies and mode shapes			16aug01
*/
void modal_results(
		FILE *fp
		, int nJ, int nM, int nI, int DoF
		, float **M, float *f, float **V
		, float total_mass, float struct_mass
		, int iter, int sumR, int modes
		, float shift, int lump, float tol
		, int ok
);

/**	
	create mesh data of deformed and undeformed mesh, use gnuplot	22feb99
	useful gnuplot options: set noxtics noytics noztics noborder view nokey
*/
void mesh(
		char IO_file[], char meshfile[], char plotfile[]
		, char *title, int nJ, int nM, int DoF
		, float *x, float *y, float *z, float *L
		, int *J1, int *J, float *p, float *D
		, float exg, int anlyz
);

/**
	create mesh data of the mode-shape meshes, use gnuplot	19oct98
	useful gnuplot options: set noxtics noytics noztics noborder view nokey
*/
void modal_mesh(
		char IO_file[], char meshfile[], char modefile[]
		, char plotfile[], char *title
		, int nJ, int nM, int DoF, int modes
		, float *x, float *y, float *z, float *L
		, int *J1, int *J2, float *p
		, float **M, float *f, float **V
		, float exg, int anlyz
);

/**
	create mesh data of animated mode-shape meshes, use gnuplot	16dec98
	useful gnuplot options: set noxtics noytics noztics noborder view nokey
	mpeg movie example:   % convert mesh_file-03-f-*.ps mode-03.mpeg
	... requires ImageMagick and mpeg2vidcodec packages
*/
void animate(
		char IO_file[], char meshfile[], char modefile[], char plotfile[]
		, char *title
		, int anim[]
		, int nJ, int nM, int DoF, int modes
		, float *x, float *y, float *z, float *L, float *p
		, int *J1, int *J2, float *f, float **V
		, float exg
		, int pan
);

/**
	computes cubic deflection functions from beam end deflections
	and beam end rotations.  Saves deflected shapes to a file.  These bent shapes
	are exact for mode-shapes, and for frames loaded at their joints.
*/
void bent_beam(
		FILE *fp, int j1, int j2
		, float *x, float *y, float *z
		, float L, float p, float *D
		, float exg
);

/** print a set of dots (periods) */
void dots(int n);

