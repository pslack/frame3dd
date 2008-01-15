#include "common.h"

#include <stdio.h>

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

/* get a line from the input file	*/
void getline (
		FILE	*fp
		, char    *s
		, int     lim
);

/* re-write input file without comments */
void parse_input(FILE *fp);

void getline_no_comment(
		FILE *fp            /**< pointer to the file from which to read */
		,char *s             /**< pointer to the string to which to write */
		,int lim            /**< the longest anticipated line length  */
);

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

void read_reactions(
		FILE *fp
		, int DoF, int *nD, int *nR
		, int nJ, float *Dp, int *R, int *sumR 
);


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

void read_condense (
		FILE *fp
		, int nJ, int modes
		, int *nC, int *Cdof, int *Cmethod, int *q, int *m
);

void control_data(
		FILE *fp
		, char *title, int nJ, int nM, int nF, int nD, int nR, int nW, int nP,int nT
		, float *x, float *y, float *z, float *r
		, int *J1, int *J2
		, float *Ax, float *Asy, float *Asz, float *J, float *Iy, float *Iz
		, float *E, float *G, float *p, float *F, float *Dp
		, float *R
		, float **W, float **P, float **T
		, int shear, int anlyz, int geom
);

void save_results (
		FILE *fp, int nJ, int nM, int DoF, int *J1, int *J2
		, float *F, float *D, int *R
		, float **Q, float err
		, int ok 
);

void modal_results(
		FILE *fp
		, int nJ, int nM, int nI, int DoF
		, float **M, float *f, float **V
		, float total_mass, float struct_mass
		, int iter, int sumR, int modes
		, float shift, int lump, float tol
		, int ok
);

void mesh(
		char IO_file[], char meshfile[], char plotfile[]
		, char *title, int nJ, int nM, int DoF
		, float *x, float *y, float *z, float *L
		, float *J1, float *J, float *p, float *D
		, float exg, int anlyz
);

void modal_mesh(
		char IO_file[], char meshfile[], char modefile[]
		, char plotfile[], char *title
		, int nJ, int nM, int DoF, int modes
		, float *x, float *y, float *z, float *L
		, float *J1, float *J2, float *p
		, float **M, float *f, float **V
		, float exg, int anlyz
);

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

void bent_beam(
		FILE *fp, int j1, int j2
		, float *x, float *y, float *z
		, float L, float p, float *D
		, float exg
);

void dots(int n);

