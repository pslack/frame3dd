/*
 This file is part of FRAME3DD:
 Static and dynamic structural analysis of 2D and 3D frames and trusses with
 elastic and geometric stiffness.
 ---------------------------------------------------------------------------
 http://frame3dd.sourceforge.net/
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
	matrix decomposition, back-substitution, and multiplication
*/

/*
 * ==========================================================================
 *
 *       Filename:  matrix.c
 *
 *    Description:  Matrix math functions
 *
 *	Version:  1.0
 *	Created:  12/30/11 18:07:41
 *       Revision:  none
 *       Compiler:  gcc
 *
 *	 Author:  Henri P. Gavin (hpgavin), h p gavin ~at~ duke ~dot~ e d v
 *	Company:  Duke Univ.
 *
 * ==========================================================================
 */


#include <math.h>
#include <stdio.h>

#include "matrix.h"
#include "nrutil.h"

/*------------------------------------------------------------------------------
LU_DCMP  -  Solves [A]{x} = {b} simply and efficiently by performing an 
  LU - decomposition of [A].  No pivoting is performed. 
  [A] is a diagonally dominant matrix of dimension [1..n][1..n]. 
  {b} is a r.h.s. vector of dimension [1..n].
  {b} is updated using [LU] and then back-substitution is done to obtain {x}.  
  {b} is replaced by {x} and [A] is replaced by the LU - reduction of itself.

  usage:  double **A, *b;
	  int   n, reduce, solve, pd; 
	  lu_dcmp ( A, n, b, reduce, solve, &pd );		     5may98

------------------------------------------------------------------------------*/
void lu_dcmp (
	double **A, /* the system matrix, and it's LU- reduction	*/
	int n,     /* the dimension of the matrix		*/	
	double *b,  /* the right hand side vector, and the solution vector  */
	int reduce,/* 1: do a forward reduction; 0: don't do the reduction */
	int solve, /* 1: do a back substitution for {x};  0: do no bk-sub'n */
	int *pd	   /* 1: positive diagonal  and  successful LU decomp'n	*/
){
	double	pivot;		/* a diagonal element of [A]		*/
	int	i, j, k;

	*pd = 1;
	if ( reduce ) {			/* forward reduction of [A]	*/

	    for (k=1; k <= n; k++) {
		if ( 0.0 == (pivot = A[k][k]) ) {
		    fprintf(stderr," lu_dcmp: zero found on the diagonal\n");
		    fprintf(stderr," A[%d][%d] = %11.4e\n", k, k, A[k][k] );
		    *pd = 0;
		    return;
		}
		for (i = k+1; i <= n; i++) {
		    A[i][k] /= pivot;
		    for (j=k+1; j <= n; j++)	A[i][j] -= A[i][k]*A[k][j];
		}
	    }
	}		/* the forward reduction of [A] is now complete	*/

	if ( solve ) {		/* back substitution to solve for {x}	*/

	    /* {b} is run through the same forward reduction as was [A]	*/

	    for (k=1; k <= n; k++)
		for (i=k+1; i <= n; i++)	b[i] -= A[i][k]*b[k];

	    /* now back substitution is conducted on {b};  [A] is preserved */

	    for (j=n; j >= 2; j--)
		for (i=1; i <= j-1; i++)	b[i] -= b[j]*A[i][j]/A[j][j];

	    /* finally we solve for the {x} vector			*/

	    for (i=1; i<=n; i++)		b[i] /= A[i][i];
	} 

	/* TAA DAAAAAAAA! {b} is now {x} and is ready to be returned	*/

	return;
}


/*----------------------------------------------------------------------------- 
LDL_DCMP  -  Solves [A]{x} = {b} simply and efficiently by performing an 
  L D L' - decomposition of [A].  No pivoting is performed.  
  [A] is a symmetric diagonally-dominant matrix of dimension [1..n][1..n]. 
  {b} is a r.h.s. vector of dimension [1..n].
  {b} is updated using L D L' and then back-substitution is done to obtain {x}. 
  {b} is returned unchanged.  ldl_dcmp(A,n,d,x,x,1,1,&pd) is valid.  
      The lower triangle of [A] is replaced by the lower triangle L of the 
      L D L' reduction.  The diagonal of D is returned in the vector {d}

 usage: double **A, *d, *b, *x;
	int   n, reduce, solve, pd;
	ldl_dcmp ( A, n, d, b, x, reduce, solve, &pd );

 H.P. Gavin, Civil Engineering, Duke University, hpgavin@duke.edu  9 Oct 2001
 Bathe, Finite Element Procecures in Engineering Analysis, Prentice Hall, 1982
-----------------------------------------------------------------------------*/
void ldl_dcmp (
	double **A,	/**< the system matrix, and L of the L D L' decomp.*/
	int n,		/**< the dimension of the matrix		*/
	double *d,	/**< diagonal of D in the  L D L' - decomp'n    */
	double *b,	/**< the right hand side vector		 */
	double *x,	/**< the solution vector			*/
	int reduce,	/**< 1: do a forward reduction of A; 0: don't   */
	int solve,	/**< 1: do a back substitution for {x}; 0: don't */
	int *pd		/**< 1: definite matrix and successful L D L' decomp'n*/
){
	int	i, j, k, m;
	*pd = 0;	/* number of negative elements on the diagonal of D */

	if ( reduce ) {		/* forward column-wise reduction of [A]	*/

	    for (j=1; j<=n; j++) {

	    	for (m=1, i=1; i < j; i++) 	/* scan the sky-line	*/
			if ( A[i][j] == 0.0 ) ++m; else	break;
		
		for (i=m; i < j; i++) {
			A[j][i] = A[i][j];
			for (k=m; k < i; k++) A[j][i] -= A[j][k]*A[i][k];
		}

		d[j] = A[j][j];
	    	for (i=m; i < j; i++)	d[j] -= A[j][i]*A[j][i]/d[i];

	    	for (i=m; i < j; i++)	A[j][i] /= d[i];
		if ( d[j] == 0.0 ) {
		    fprintf(stderr," ldl_dcmp(): zero found on diagonal ...\n");
		    fprintf(stderr," d[%d] = %11.4e\n", j, d[j] );
		    return;
		}
		if ( d[j] < 0.0 ) (*pd)--;
	    }
	}		/* the forward reduction of [A] is now complete	*/

	if ( solve ) {		/* back substitution to solve for {x}   */

		/* {x} is run through the same forward reduction as was [A] */

	    for (i=1; i <= n; i++) {
		x[i] = b[i];
		for (j=1; j < i; j++)		x[i] -= A[i][j]*x[j];
	    }

	    for (i=1; i <= n; i++)		x[i] /= d[i];

	    /* now back substitution is conducted on {x};  [A] is preserved */

	    for (i=n; i > 1; i--) 
		for (j=1; j < i; j++)		x[j] -= A[i][j]*x[i];

	} 
	return;
}


/*----------------------------------------------------------------------------
LDL_MPROVE  -  Improves a solution vector x[1..n] of the linear set of equations
 [A]{x} = {b}.  The matrix A[1..n][1..n], and the vectors b[1..n] and x[1..n]
 are input, as is the dimension n.   The matrix [A] is the L D L'
 decomposition of the original system matrix, as returned by ldl_dcmp().
 Also input is the diagonal vector, {d} of [D] of the L D L' decompositon.
 On output, only {x} is modified to an improved set of values.

 usage: double **A, *d, *b, *x, err;
	int   n, ok;
	ldl_mprove ( A, n, d, b, x, &err, &ok );

 H.P. Gavin, Civil Engineering, Duke University, hpgavin@duke.edu  4 May 2001
-----------------------------------------------------------------------------*/
void ldl_mprove(
	double **A, int n, double *d, double *b, double *x, 
	double *err, int *ok
){
	double  sdp;		/* accumulate the r.h.s. in double precision */
	double   *r,		/* the residual error		  	*/
		err_new=0,	/* the RMS error of the solution	*/
		*dvector();	/* allocate memory for a vector	of doubles */
	int	j,i, pd;
	void	ldl_dcmp(),
		free_dvector();

	r=dvector(1,n);

	for (i=1;i<=n;i++) {		/* calculate the r.h.s. of      */
		sdp = b[i];		/* [A]{r} = {b} - [A]{x+r}      */
		for (j=1;j<=n;j++) {	/* A in upper triangle only     */
			if ( i <= j )   sdp -= A[i][j] * x[j];
			else		sdp -= A[j][i] * x[j];
		}
		r[i] = sdp;
	}

	ldl_dcmp ( A, n, d, r, r, 0, 1, &pd );   /* solve for the error term */

	for (i=1;i<=n;i++)	err_new += r[i]*r[i];

	err_new = sqrt ( err_new / (double) n );

	*ok = 0;
	if ( err_new / *err < 0.90 ) {		/* good improvement */
				/* subtract the error from the old solution */
		for (i=1;i<=n;i++)	x[i] += r[i];
		*err = err_new;
		*ok = 1;	/* the solution has improved		*/
	}

	free_dvector(r,1,n);
	return;
}


/*----------------------------------------------------------------------------
PSEUDO_INV - calculate the pseudo-inverse of A ,
	     Ai = inv ( A'*A + beta * trace(A'*A) * I ) * A' 
	     beta is a regularization factor, which should be small (1e-10)
	     A is m by n      Ai is m by n			      8oct01
-----------------------------------------------------------------------------*/
void pseudo_inv(
	double **A, double **Ai, int n, int m, double beta, int verbose
){
	double	*diag, *b, *x, **AtA, **AtAi, tmp, tr_AtA=0.0,
		*dvector(), **dmatrix(), error;
	int     i,j,k, ok; 
	void	ldl_dcmp(), ldl_mprove(), free_dvector(), free_dmatrix();

	diag = dvector(1,n);
	b    = dvector(1,n);
	x    = dvector(1,n);
	AtA  = dmatrix(1,n,1,n);
	AtAi = dmatrix(1,n,1,n);

	if (beta>1) fprintf(stderr," pseudo_inv: warning beta = %lf\n", beta);

	for (i=1; i<=n; i++) {
		diag[i] = x[i] = b[i] = 0.0;
		for (j=i; j<=n; j++)    AtA[i][j] = AtA[j][i] = 0.0;
	}

	for (i=1; i<=n; i++) {			/* compute A' * A */
		for (j=1; j<=n; j++) {
			tmp = 0.0;
			for (k=1; k<=m; k++) tmp += A[k][i] * A[k][j];
			AtA[i][j] = tmp;
		}
	}
	for (i=1; i<=n; i++)			    /* make symmetric */
		for (j=i; j<=n; j++)
			AtA[i][j]=AtA[j][i] = 0.5*(AtA[i][j] + AtA[j][i]);

	tr_AtA = 0.0;
	for (i=1; i<=n; i++) tr_AtA += AtA[i][i];       /* trace of AtA */
	for (i=1; i<=n; i++) AtA[i][i] += beta*tr_AtA;	/* add beta I */

	ldl_dcmp ( AtA, n, diag, b, x, 1, 0, &ok );	/*  L D L'  decomp */

	for (j=1; j<=n; j++) {				/* compute inv(AtA) */

		for (k=1; k<=n; k++)  b[k] = 0.0;
		b[j] = 1.0;
		ldl_dcmp( AtA, n, diag, b, x, 0, 1, &ok ); /* L D L' bksbtn */

		if ( verbose )
		 fprintf(stdout,"  RMS matrix error:"); /*improve the solution*/
		error = 1.0; ok = 1;
		do {
			ldl_mprove ( AtA, n, diag, b, x, &error, &ok );
			if ( verbose ) fprintf(stdout,"%9.2e", error );
		} while ( ok );
		if ( verbose ) fprintf(stdout,"\n");

		for (k=1; k<=n; k++)  AtAi[k][j] = x[k];  /* save inv(AtA) */
	}

	for (i=1; i<=n; i++)			    /* make symmetric */
		for (j=i; j<=n; j++)
			AtAi[i][j]=AtAi[j][i] = 0.5*(AtAi[i][j] + AtAi[j][i]);

	for (i=1; i<=n; i++) {			/* compute inv(A'*A)*A'	*/
		for (j=1; j<=m; j++) {
			tmp = 0.0;
			for (k=1; k<=n; k++)    tmp += AtAi[i][k]*A[j][k];
			Ai[i][j] = tmp;
		}
	}

	free_dmatrix (AtAi, 1,n,1,n);
	free_dmatrix (AtA,  1,n,1,n);
	free_dvector (x, 1,n);
	free_dvector (b, 1,n);
	free_dvector (diag, 1,n);

	return;
}


/* ----------------------------------------------------------------------------
 * PRODABj -  matrix-matrix multiplication for symmetric A	      27apr01
 *		 u = A * B(:,j)
 * --------------------------------------------------------------------------*/
void prodABj ( double **A, double **B, double *u, int n, int j )
{
	int     i, k;

	for (i=1; i<=n; i++)    u[i] = 0.0;

	for (i=1; i<=n; i++) {
		for (k=1; k<=n; k++) {
			if ( i <= k )   u[i] += A[i][k]*B[k][j];
			else	    u[i] += A[k][i]*B[k][j];
		}
	}
	return;
}


/* ----------------------------------------------------------------------------
 * prodAB - matrix-matrix multiplication      C = A * B			27apr01
 * --------------------------------------------------------------------------*/
void prodAB ( double **A, double **B, double **C, int I, int J, int K )
{
	int     i, j, k;

	for (i=1; i<=I; i++)
		for (k=1; k<=K; k++)
			C[i][k] = 0.0;

	for (i=1; i<=I; i++)
		for (k=1; k<=K; k++)
			for (j=1; j<=J; j++)
				C[i][k] += A[i][j]*B[j][k];
	return;
}


/* ---------------------------------------------------------------------------
INVAB  -  calculate product inv(A) * B  
	 A is n by n      B is n by m				    6jun07
----------------------------------------------------------------------------*/
void invAB( double **A, double **B, int n, int m, double **AiB, int *ok, int verbose)
{
	double  *diag, *b, *x, error;
	int     i,j,k;

	diag = dvector(1,n);
	x    = dvector(1,n);
	b    = dvector(1,n);

	for (i=1; i<=n; i++) diag[i] = x[i] = 0.0;

	ldl_dcmp( A, n, diag, b, x, 1, 0, ok );	 /*   L D L'  decomp */
	if ( *ok < 0 ) {
		fprintf(stderr," Make sure that all six");
		fprintf(stderr," rigid body translations are restrained!\n");
	}

	for (j=1; j<=m; j++) {

		for (k=1; k<=n; k++)  b[k] = B[k][j];
		ldl_dcmp( A, n, diag, b, x, 0, 1, ok ); /*   L D L'  bksbtn */

		if ( verbose ) fprintf(stdout,"    LDL' RMS matrix precision:");
		error = *ok = 1;
		do {				    /* improve the solution*/
			ldl_mprove ( A, n, diag, b, x, &error, ok );
			if ( verbose ) fprintf(stdout,"%9.2e", error );
		} while ( *ok );
		if ( verbose ) fprintf(stdout,"\n");

		for (i=1; i<=n; i++)    AiB[i][j] = x[i];
	}

	free_dvector(diag,1,n);
	free_dvector(x,1,n);
	free_dvector(b,1,n);
}


/* ---------------------------------------------------------------------------
 * XTinvAY  -  calculate quadratic form with inverse matrix   X' * inv(A) * Y   
 *	   A is n by n    X is n by m     Y is n by m		    15sep01
 * --------------------------------------------------------------------------*/
void xtinvAy(
	double **X, double **A, double **Y, int n, int m, double **Ac, int verbose
){
	double  *diag, *x, *y, error;
	int     i,j,k, ok;

	diag = dvector(1,n);
	x    = dvector(1,n);
	y    = dvector(1,n);

	for (i=1; i<=n; i++) diag[i] = x[i] = 0.0;

	ldl_dcmp( A, n, diag, y, x, 1, 0, &ok );	/*   L D L'  decomp */

	for (j=1; j<=m; j++) {

		for (k=1; k<=n; k++)  y[k] = Y[k][j];
		ldl_dcmp( A, n, diag, y, x, 0, 1, &ok ); /*   L D L'  bksbtn */

		if ( verbose ) fprintf(stdout,"    LDL' RMS matrix precision:");
		error = ok = 1;
		do {				    /* improve the solution*/
			ldl_mprove ( A, n, diag, y, x, &error, &ok );
			if ( verbose ) fprintf(stdout,"%9.2e", error );
		} while ( ok );
		if ( verbose ) fprintf(stdout,"\n");

		for (i=1; i<=m; i++) {
			Ac[i][j] = 0.0;
			for (k=1; k<=n; k++) Ac[i][j] += X[k][i] * x[k];
		}
	}

	free_dvector(diag,1,n);
	free_dvector(x,1,n);
	free_dvector(y,1,n);
}


/*-----------------------------------------------------------------------------
 * xtAx - carry out matrix-matrix-matrix multiplication for symmetric A  7nov02
 *	 C = X' A X     C is J by J      X is N by J     A is N by N      
 * ---------------------------------------------------------------------------*/
void xtAx(double **A, double **X, double **C, int N, int J)
{

	double  **AX;
	int     i,j,k;

	AX = dmatrix(1,N,1,J);

	for (i=1; i<=J; i++)    for (j=1; j<=J; j++)     C[i][j] = 0.0;
	for (i=1; i<=N; i++)    for (j=1; j<=J; j++)    AX[i][j] = 0.0;
		
	for (i=1; i<=N; i++) {	  /*  use upper triangle of A */
		for (j=1; j<=J; j++) {
			for (k=1; k<=N; k++) {
				if ( i <= k )   AX[i][j] += A[i][k] * X[k][j];
				else	    AX[i][j] += A[k][i] * X[k][j];
			}
		} 
	}       
		
	for (i=1; i<=J; i++)    
		for (j=1; j<=J; j++)
			for (k=1; k<=N; k++) 
				C[i][j] += X[k][i] * AX[k][j];
			
	for (i=1; i<=J; i++)	    /*  make  C  symmetric */
		for (j=i; j<=J; j++)
			C[i][j] = C[j][i] = 0.5 * ( C[i][j] + C[j][i] );

	free_dmatrix(AX,1,N,1,J);
	return;
}


/* ----------------------------------------------------------------------------
 * xtAy - carry out vector-matrix-vector multiplication for symmetric A  7apr94
 * ------------------------------------------------------------------------- */
double xtAy(double *x, double **A, double *y, int n, double *d)
{

	double  xtAy = 0.0;
	int     i,j;

	for (i=1; i<=n; i++) {				/*  d = A y  */
		d[i]  = 0.0;    
		for (j=1; j<=n; j++) {	  //  A in upper triangle only 
			if ( i <= j )   d[i] += A[i][j] * y[j];
			else	    d[i] += A[j][i] * y[j];
		}
	}
	for (i=1; i<=n; i++)    xtAy += x[i] * d[i];	/*  xAy = x' A y  */
	return ( xtAy );
}


/* ---------------------------------------------------------------------------
 * invAXinvA -  calculate quadratic form with inverse matrix 
 *	   replace X with inv(A) * X * inv(A) 
 *	   A is n by n and symmetric   X is n by n and symmetric    15sep01
 * -------------------------------------------------------------------------*/
void invAXinvA ( double **A, double **X, int n, int verbose )
{
	double  *diag, *b, *x, **Ai, **XAi, Aij, error;
	int     i,j,k, ok;

	diag = dvector(1,n);
	x    = dvector(1,n);
	b    = dvector(1,n);
	Ai   = dmatrix(1,n,1,n);
	XAi  = dmatrix(1,n,1,n);

	for (i=1; i<=n; i++) {
		diag[i] = x[i] = b[i] = 0.0;
		for (j=1; j<=n; j++)    XAi[i][j] = Ai[i][j] = 0.0;
	}

	ldl_dcmp ( A, n, diag, b, x, 1, 0, &ok );       /*   L D L'  decomp */

	for (j=1; j<=n; j++) {			  /*  compute inv(A)  */

		for (k=1; k<=n; k++)  b[k] = 0.0;
		b[j] = 1.0;
		ldl_dcmp( A, n, diag, b, x, 0, 1, &ok ); /*   L D L'  bksbtn */

		if ( verbose ) fprintf(stdout,"    LDL' RMS matrix precision:");
		error = ok = 1;
		do {				    /* improve the solution*/
			ldl_mprove ( A, n, diag, b, x, &error, &ok );
			if ( verbose ) fprintf(stdout,"%9.2e", error );
		} while ( ok );
		if ( verbose ) fprintf(stdout,"\n");

		for (k=1; k<=n; k++)  Ai[j][k] = x[k];  /*  save inv(A) */
	}

	for (i=1; i<=n; i++)			    /*  make symmetric */
		for (j=i; j<=n; j++)
			Ai[i][j] = Ai[j][i] = 0.5 * ( Ai[i][j] + Ai[j][i] );

	for (i=1; i<=n; i++) {		  /*  compute X * inv(A)   */
		for (j=1; j<=n; j++) {
			Aij = 0.0;
			for (k=1; k<=n; k++)    Aij += X[i][k]*Ai[k][j];
			XAi[i][j] = Aij;
		}
	}

	for (i=1; i<=n; i++) {	  /*  compute inv(A) * X * inv(A)  */
		for (j=1; j<=n; j++) {
			Aij = 0.0;
			for (k=1; k<=n; k++)    Aij += Ai[i][k] * XAi[k][j];
			X[i][j] = Aij;
		}
	}
	for (i=1; i<=n; i++)			    /*  make symmetric */
		for (j=i; j<=n; j++)
			X[i][j] = X[j][i] = 0.5 * ( X[i][j] + X[j][i] );

	free_dvector ( diag, 1,n );
	free_dvector ( x, 1,n );
	free_dvector ( b, 1,n );
	free_dmatrix ( Ai, 1,n,1,n );
	free_dmatrix ( XAi, 1,n,1,n );

}


/* ---------------------------------------------------------------------------
 *  REL_NORM -  compute the relative 2-norm between two vectors       26dec01 
 *       compute the relative 2-norm between two vectors N and D
 *	       return  ( sqrt(sum(N[i]*N[i]) / sqrt(D[i]*D[i]) )
 *
 * ------------------------------------------------------------------------- */
double rel_norm( double *N, double *D, int n )
{
	double  nN = 0.0, nD = 0.0;
	int     i;

	for (i=1; i<=n; i++)    nN += ( N[i]*N[i] );
	for (i=1; i<=n; i++)    nD += ( D[i]*D[i] );

	return ( sqrt(nN) / sqrt(nD) );
}


