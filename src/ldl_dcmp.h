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
 LDL_DCMP  
 Solves [A]{x} = {b} simply and efficiently by performing an 
 L D L' - decomposition of [A].  No pivoting is performed.  
 [A] is a symmetric diagonally-dominant matrix of dimension [1..n][1..n].
 {b} is a r.h.s. vector of dimension [1..n].
 {b} is updated using L D L' and then back-substitution is done to obtain {x}. 
 {b} is returned unchanged.  ldl_dcmp(A,n,d,x,x,1,1,&pd) is valid.  
     The lower triangle of [A] is replaced by the lower triangle L of the 
     L D L' reduction.  The diagonal of D is returned in the vector {d}
*/
void ldl_dcmp (
	double **A,	/**< the system matrix, and L of the L D L' decomp.*/
	int n,		/**< the dimension of the matrix		*/	
	double *d,	/**< diagonal of D in the  L D L' - decomp'n	*/
	double *b,	/**< the right hand side vector			*/
	double *x,	/**< the solution vector			*/
	int reduce,	/**< 1: do a forward reduction of A; 0: don't	*/
	int solve,	/**< 1: do a back substitution for {x}; 0: don't */
	int *pd		/**< 1: definite matrix and successful L D L' decomp'n*/
);

/**
	Improves a solution vector x[1..n] of the linear set of equations
	[A]{x} = {b}. 
	The matrix A[1..n][1..n], and the vectors b[1..n] and x[1..n]
	are input, as is the dimension n.   The matrix [A] is the L D L'
	decomposition of the original system matrix, as returned by ldl_dcmp().
	Also input is the diagonal vector, {d} of [D] of the L D L' decompositon.
	On output, only {x} is modified to an improved set of values.
*/
void ldl_mprove(
	double **A, 
	int n, double *d, double *b, double *x,
	double *err,	/**< the RMS error of the solution		*/
	int *ok
);

/**
	calculate the pseudo-inverse of A ,
	Ai = inv ( A'*A + beta * trace(A'*A) * I ) * A' 
       	beta is a regularization factor, which should be small (1e-10)
	A is m by n      Ai is m by n  
*/
void pseudo_inv(
	double **A,	/**< an n-by-m matrix				*/
	double **Ai,	/**< the pseudo-inverse of A			*/
	int n, int m,	/**< matrix dimensions				*/
	double beta	/**< regularization factor			*/
);

#endif /* FRAME_LDL_DCMP_H */

