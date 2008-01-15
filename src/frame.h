#ifndef FRAME_FRAME_H
#define FRAME_FRAME_H

/* note, common.h includes a sneaky #define that rewrites 'float' to 'double'. */
#include "common.h"

/* form the global stiffness matrix	*/
void assemble_K(
	float **K
	, int DoF, int nM
	, float *x, float *y, float *z, float *r, float *L, float *Le
	, int *J1, int *J2
	, float *Ax, float *Asy, float *Asz
	, float *J, float *Iy, float *Iz, float *E, float *G, float *p
	, int shear, int geom, float **Q
);


void atma();		/* carry out the coordinate transfm'n	*/

void
	apply_reactions(), /* apply boundary conditions		*/
	solve_system(),	/* solve a linear system via LDL' dcmp	*/
	end_forces(),	/* evaluate the member end forces	*/
	equilibrium(),	/* perform an equilibrium check		*/
	assemble_M(),	/* form the global mass matrix		*/
	condense(),	/* static matrix condensation		*/
	guyan(),	/* Guyan reduction of matrices Md , Kd	*/
	dyn_conden(),	/* dynamic condensation of Md and Kd	*/
	stodola(),	/* lower generalized eigenval & eigenvec*/
	subspace(),	/* lower generalized eigenval & eigenvec*/
	dots(),		/* print a string of dots (periods)	*/
	deallocate();	/* release the allocated memory		*/

#endif /* FRAME_FRAME_H */

