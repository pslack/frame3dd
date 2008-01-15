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

