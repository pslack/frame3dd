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
	LDL' decomposition and back-substitution
*/
#ifndef FRAME_LDL_DCMP_H
#define FRAME_LDL_DCMP_H

/**
	L D L' decompositon and back-substitution 
*/
void ldl_dcmp (
		double **A /**< the system matrix, and L of the L D L' decomposition	*/
		, int n /**< the dimension of the matrix */	
		, double *d /**< diagonal of D in the  L D L' - decomposition */
		, double *b /**< the right hand side vector */
		, double *x /**< the solution vector */
		, int reduce /**< 1: do a forward reduction of A; 0: do no reduction */
		, int solve /**< 1: do a back substitution for {x}; 0: do no bk-sub'n */
		, int *pd /**< 1: definite matrix  and  successful L D L' decomp'n	*/
);

/**
	iterative improvement to the solution
*/
void ldl_mprove(
		 double **A
		, int n, double *d, double *b, double *x
		, double *err, int *ok
);

/**
	pseudo-inverse
*/
void pseudo_inv(
		double **A, double **Ai
		, int n, int m, double beta
);

/**
	relative norm
*/
double rel_norm(double *N, double *D, int n);

#endif /* FRAME_LDL_DCMP_H */

