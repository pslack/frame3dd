/* 
 ---------------------------------------------------------------------------
 FILE   ldl_dcmp.c - routines to perform L D L' - decomposition
 ---------------------------------------------------------------------------
 Copyright (C) 1992-2002  Henri P. Gavin
 
 This program is free software; you may redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 http://www.fsf.org/copyleft/gpl.html
 
 You should have received a copy of the GNU General Public License, gpl.txt,
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 ---------------------------------------------------------------------------
 Henri P. Gavin                                             hpgavin@duke.edu
 Department of Civil and Environmental Engineering
 Duke University, Box 90287
 Durham, NC  27708--0287
 ---------------------------------------------------------------------------
*/


#include <stdio.h>
#include <math.h>

/* ----------------------- double precision -------------------------------- */
#define float double
#define vector dvector
#define matrix dmatrix
#define free_vector free_dvector
#define free_matrix free_dmatrix
#define save_matrix save_dmatrix
#define show_matrix show_dmatrix
#define show_vector show_dvector
/* ----------------------- double precision -------------------------------- */

/*----------------------------------------------------------------------------- 
LDL_DCMP  -  Solves [A]{x} = {b} simply and efficiently by performing an 
  L D L' - decomposition of [A].  No pivoting is performed.  
  [A] is a symmetric diagonally-dominant matrix of dimension [1..n][1..n]. 
  {b} is a r.h.s. vector of dimension [1..n].
  {b} is updated using L D L' and then back-substitution is done to obtain {x}. 
  {b} is returned unchanged.  ldl_dcmp(A,n,d,x,x,1,1,&pd) is valid.  
      The lower triangle of [A] is replaced by the lower triangle L of the 
      L D L' reduction.  The diagonal of D is returned in the vector {d}

 usage: float **A, *d, *b, *x;
	int   n, reduce, solve, pd;
	ldl_dcmp ( A, n, d, b, x, reduce, solve, &pd );

 H.P. Gavin, Civil Engineering, Duke University, hpgavin@duke.edu  9 Oct 2001
 Bathe, Finite Element Procecures in Engineering Analysis, Prentice Hall, 1982
-----------------------------------------------------------------------------*/
void ldl_dcmp ( A, n, d, b, x, reduce, solve, pd )
float	**A,	/* the system matrix, and L of the L D L' decomposition	*/
	*d,	/* diagonal of D in the  L D L' - decomposition		*/
	*b,	/* the right hand side vector				*/
	*x;	/* the solution vector					*/
int	n,	/* the dimension of the matrix				*/	
	reduce,	/* 1: do a forward reduction of A;    0: do no reduction */
	solve,	/* 1: do a back substitution for {x}; 0: do no bk-sub'n	*/
	*pd;	/* 1: definite matrix  and  successful L D L' decomp'n	*/
{
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

 usage: float **A, *d, *b, *x, err;
	int   n, ok;
	ldl_mprove ( A, n, d, b, x, &err, &ok );

 H.P. Gavin, Civil Engineering, Duke University, hpgavin@duke.edu  4 May 2001
-----------------------------------------------------------------------------*/
void ldl_mprove ( A, n, d, b, x, err, ok )
float   **A, *d, *b, *x, *err;
int     n, *ok;
{
	double  sdp;		/* accumulate the r.h.s. in double precision */
	float   *r,		/* the residual error		  	*/
		err_new=0,	/* the RMS error of the solution	*/
		*vector();	/* allocate memory for a vector	of floats */
	int	j,i, pd;
	void	ldl_dcmp(),
		free_vector();

	r=vector(1,n);

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

	err_new = sqrt ( err_new / (float) n );

	*ok = 0;
	if ( err_new / *err < 0.90 ) {		/* good improvement */
				/* subtract the error from the old solution */
		for (i=1;i<=n;i++)	x[i] += r[i];
		*err = err_new;
		*ok = 1;	/* the solution has improved		*/
	}

	free_vector(r,1,n);
	return;
}


/*----------------------------------------------------------------------------
PSEUDO_INV - calculate the pseudo-inverse of A ,
             Ai = inv ( A'*A + beta * trace(A'*A) * I ) * A' 
             beta is a regularization factor, which should be small (1e-10)
             A is m by n      Ai is m by n                              8oct01
-----------------------------------------------------------------------------*/
void pseudo_inv ( A, Ai, n, m, beta )
float   **A, **Ai, beta; 
int     n, m;
{
	float	*diag, *b, *x, **AtA, **AtAi, tmp, tr_AtA=0.0,
		*vector(), **matrix(), error;
	int     i,j,k, ok, disp=0;
	void	ldl_dcmp(), ldl_mprove(), free_vector(), free_matrix();

        diag = vector(1,n);
        b    = vector(1,n);
        x    = vector(1,n);
        AtA  = matrix(1,n,1,n);
	AtAi = matrix(1,n,1,n);

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
	for (i=1; i<=n; i++)                            /* make symmetric */
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

		if (disp)
		 fprintf(stderr,"  RMS matrix error:"); /*improve the solution*/
		error = 1.0; ok = 1;
		do {
			ldl_mprove ( AtA, n, diag, b, x, &error, &ok );
			if (disp) fprintf(stderr,"%9.2e", error );
		} while ( ok );
		if (disp) fprintf(stderr,"\n");

		for (k=1; k<=n; k++)  AtAi[k][j] = x[k];  /* save inv(AtA) */
        }

	for (i=1; i<=n; i++)                            /* make symmetric */
		for (j=i; j<=n; j++)
			AtAi[i][j]=AtAi[j][i] = 0.5*(AtAi[i][j] + AtAi[j][i]);

	for (i=1; i<=n; i++) {			/* compute inv(A'*A)*A'	*/
		for (j=1; j<=m; j++) {
			tmp = 0.0;
			for (k=1; k<=n; k++)    tmp += AtAi[i][k]*A[j][k];
			Ai[i][j] = tmp;
		}
	}

	free_matrix (AtAi, 1,n,1,n);
	free_matrix (AtA,  1,n,1,n);
	free_vector (x, 1,n);
	free_vector (b, 1,n);
	free_vector (diag, 1,n);

	return;
}


/* ---------------------------------------------------------------------------
REL_NORM -  compute the relative 2-norm between two vectors       26dec01 
--------------------------------------------------------------------------- */
float rel_norm( N, D, n )
float	*N, *D;
int	n;
{
	float	nN = 0.0, nD = 0.0;
	int	i;

	for (i=1; i<=n; i++)	nN += N[i]*N[i];
	for (i=1; i<=n; i++)	nD += D[i]*D[i];

	return  ( sqrt(nN) / sqrt(nD) );
}
