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
*/
#ifndef FRAME_LDL_DCMP_H
#define FRAME_LDL_DCMP_H

#include "common.h"

/**
	L D L' decompositon and back-substitution 
*/
void ldl_dcmp (
		float **A /**< the system matrix, and L of the L D L' decomposition	*/
		, int n /**< the dimension of the matrix */	
		, float *d /**< diagonal of D in the  L D L' - decomposition */
		, float *b /**< the right hand side vector */
		, float *x /**< the solution vector */
		, int reduce /**< 1: do a forward reduction of A; 0: do no reduction */
		, int solve /**< 1: do a back substitution for {x}; 0: do no bk-sub'n */
		, int *pd /**< 1: definite matrix  and  successful L D L' decomp'n	*/
);

/**
	iterative improvement to the solution
*/
void ldl_mprove(
		 float **A
		, int n, float *d, float *b, float *x
		, float *err, int *ok
);

void pseudo_inv(
		float **A, float **Ai
		, int n, int m, float beta
);

float rel_norm(float *N, float *D, int n);

#endif /* FRAME_LDL_DCMP_H */

