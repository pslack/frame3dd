/* 
 ---------------------------------------------------------------------------
 FILE lu_dcmp.c - routins to perform L U - decomposition 
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


/*------------------------------------------------------------------------------
LU_DCMP  -  Solves [A]{x} = {b} simply and efficiently by performing an 
  LU - decomposition of [A].  No pivoting is performed. 
  [A] is a diagonally dominant matrix of dimension [1..n][1..n]. 
  {b} is a r.h.s. vector of dimension [1..n].
  {b} is updated using [LU] and then back-substitution is done to obtain {x}.  
  {b} is replaced by {x} and [A] is replaced by the LU - reduction of itself.

  usage:  float **A, *b;
	  int   n, reduce, solve, pd; 
	  lu_dcmp ( A, n, b, reduce, solve, &pd );                     5may98

------------------------------------------------------------------------------*/
void lu_dcmp ( A, n, b, reduce, solve, pd )
float	**A,	/* the system matrix, and it's LU- reduction		*/
	*b;	/* the right hand side vector, and the solution vector	*/
int	n,	/* the dimension of the matrix				*/	
	reduce,	/* 1: do a forward reduction; 0: don't do the reduction */
	solve,	/* 1: do a back substitution for {x};  0: do no bk-sub'n */
	*pd;	/* 1: positive diagonal  and  successful LU decomp'n	*/
{
	float	pivot;		/* a diagonal element of [A]		*/
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
