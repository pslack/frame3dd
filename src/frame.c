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

/* NOTE main 'driver' routine is now moved to main.c */

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "nrutil.h"

/* must come after the above, because of the sneaky #defines in common.h */
#include "frame.h"
#include "common.h"
#include "coordtrans.h"
#include "ldl_dcmp.h"

/* forward declarations */

static void elastic_K(
	float **k, float *x, float *y, float *z, float *r
	, float L, float Le
	, int j1, int j2
	, float Ax, float Asy, float Asz
	, float J, float Iy, float Iz, float E, float G, float p
	, int shear
);

static void geometric_K(
	float **k, float *x, float *y, float *z, float *r
	, float L, float Le
	, int j1, int j2, float Ax, float Asy, float Asz, float J
	, float Iy, float Iz, float E, float G, float p, float T
	, int shear
);

static void member_force(
	float *s, int M, float *x, float *y, float *z, float L, float Le
	, int j1, int j2, float Ax, float Asy, float Asz, float J
	, float Iy, float Iz, float E, float G, float p, float *D
	, int shear, int geom
);

static void lumped_M(
	float **m, float *x, float *y, float *z
	, float L, int j1, int j2, float Ax, float J, float Iy, float Iz
	, float d, float p, float BMs
);

static void consistent_M(
	float **m, float *x, float *y, float *z, float *r, float L
	, int j1, int j2, float Ax, float J, float Iy, float Iz, float d
	, float BMs, float p
);

static void invAB(
	float **A, float **B
	, int n, int m, float **AiB
	, int *ok
);

static void xtaiy(
	float **X, float **A, float **Y
	, int n, int m, float **Ac
);

/*------------------------------------------------------------------------------
ASSEMBLE_K  -  assemble global stiffness matrix from individual elements 23feb94
------------------------------------------------------------------------------*/
void assemble_K(
	float **K
	, int DoF, int nM
	, float *x, float *y, float *z, float *r, float *L, float *Le
	, int *J1, int *J2
	, float *Ax, float *Asy, float *Asz
	, float *J, float *Iy, float *Iz, float *E, float *G, float *p
	, int shear, int geom, float **Q
){
	float	**k;		/* element stiffness matrix in global coord */
	int	**ind,		/* member-structure DoF index table	*/
		i, j, ii, jj, l, ll;

	for (i=1; i<=DoF; i++)	for (j=1; j<=DoF; j++)	K[i][j] = 0.0;

	k   =  matrix(1,12,1,12);
	ind = imatrix(1,12,1,nM);


	for ( i=1; i<= nM; i++ ) {
		ind[1][i] = 6*J1[i] - 5;	ind[7][i]  = 6*J2[i] - 5;
		ind[2][i] = ind[1][i] + 1;	ind[8][i]  = ind[7][i] + 1;
		ind[3][i] = ind[1][i] + 2;	ind[9][i]  = ind[7][i] + 2;
		ind[4][i] = ind[1][i] + 3;	ind[10][i] = ind[7][i] + 3;
		ind[5][i] = ind[1][i] + 4;	ind[11][i] = ind[7][i] + 4;
		ind[6][i] = ind[1][i] + 5;	ind[12][i] = ind[7][i] + 5;
	}

	for ( i = 1; i <= nM; i++ ) {

		elastic_K ( k, x,y,z, r, L[i], Le[i], J1[i], J2[i],
		Ax[i],Asy[i],Asz[i], J[i], Iy[i],Iz[i], E[i],G[i], p[i], shear);

		if (geom)
		 geometric_K( k, x,y,z,r, L[i], Le[i], J1[i], J2[i],
		           Ax[i], Asy[i],Asz[i], 
                           J[i], Iy[i], Iz[i], 
                           E[i],G[i], p[i], -Q[i][1], shear);

		for ( l=1; l <= 12; l++ ) {
			ii = ind[l][i];
			for ( ll=1; ll <= 12; ll++ ) {
				jj = ind[ll][i];
				K[ii][jj] += k[l][ll];
			}
		}
	}
	free_matrix ( k,1,12,1,12);
	free_imatrix(ind,1,12,1,nM);
	return;
}


/*------------------------------------------------------------------------------
ELASTIC_K - space frame elastic stiffness matrix in global coordnates	22oct02
------------------------------------------------------------------------------*/
void elastic_K(
	float **k, float *x, float *y, float *z, float *r
	, float L, float Le
	, int j1, int j2
	, float Ax, float Asy, float Asz
	, float J, float Iy, float Iz, float E, float G, float p
	, int shear
){
	float   t1, t2, t3, t4, t5, t6, t7, t8, t9,     /* coord Xformn */
		Ksy, Ksz;		/* shear deformatn coefficients	*/
	int     i, j;

	coord_trans ( x, y, z, L, j1, j2,
				&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p );

	for (i=1;i<=12;i++)	for (j=1;j<=12;j++)	k[i][j] = 0.0;

	if ( shear ) {
		Ksy = 12.*E*Iz / (G*Asy*Le*Le);
		Ksz = 12.*E*Iy / (G*Asz*Le*Le);
	} else	Ksy = Ksz = 0.0;

	k[1][1]  = k[7][7]   = E*Ax / Le;
	k[2][2]  = k[8][8]   = 12.*E*Iz / ( Le*Le*Le*(1.+Ksy) );
	k[3][3]  = k[9][9]   = 12.*E*Iy / ( Le*Le*Le*(1.+Ksz) );
	k[4][4]  = k[10][10] = G*J / Le;
	k[5][5]  = k[11][11] = (4.+Ksz)*E*Iy / ( Le*(1.+Ksz) );
	k[6][6]  = k[12][12] = (4.+Ksy)*E*Iz / ( Le*(1.+Ksy) );

	k[5][3]  = k[3][5]   = -6.*E*Iy / ( Le*Le*(1.+Ksz) );
	k[6][2]  = k[2][6]   =  6.*E*Iz / ( Le*Le*(1.+Ksy) );
	k[7][1]  = k[1][7]   = -k[1][1];

	k[12][8] = k[8][12]  =  k[8][6] = k[6][8] = -k[6][2];
	k[11][9] = k[9][11]  =  k[9][5] = k[5][9] = -k[5][3];
	k[10][4] = k[4][10]  = -k[4][4];
	k[11][3] = k[3][11]  =  k[5][3];
	k[12][2] = k[2][12]  =  k[6][2];

	k[8][2]  = k[2][8]   = -k[2][2];
	k[9][3]  = k[3][9]   = -k[3][3];
	k[11][5] = k[5][11]  = (2.-Ksz)*E*Iy / ( Le*(1.+Ksz) );
	k[12][6] = k[6][12]  = (2.-Ksy)*E*Iz / ( Le*(1.+Ksy) );

/*	save_matrix ( 12, 12, k, "ke" ); /* element elastic stiffness matrix */

	atma ( t1,t2,t3,t4,t5,t6,t7,t8,t9, k, r[j1],r[j2] );	/* globalize */

	/* check and enforce symmetry */ 


	for (i=1; i<=12; i++)
	    for (j=i+1; j<=12; j++) 
		if ( k[i][j] != k[j][i] ) {
		    if (fabs(k[i][j]/k[j][i]-1.0) > 1.0e-6 &&
		       (fabs(k[i][j]/k[i][i]) > 1e-6 || 
                        fabs(k[j][i]/k[i][i]) > 1e-6) ) {
		     fprintf(stderr,"elastic_K: element stiffness matrix not symetric ...\n" ); 
		     fprintf(stderr," ... k[%d][%d] = %15.6e \n",i,j,k[i][j] ); 
		     fprintf(stderr," ... k[%d][%d] = %15.6e   ",j,i,k[j][i] ); 
		     fprintf(stderr," ... relative error = %e \n",  fabs(k[i][j]/k[j][i]-1.0) ); 
		     fprintf(stderr," ... element matrix saved in file 'kt'\n");
		     save_matrix ( 12, 12, k, "kt" ); 
		    }
		    k[i][j] = k[j][i] = 0.5 * ( k[i][j] + k[j][i] );
		}


/*	save_matrix ( 12, 12, k, "ket" );   /* transformed element matx */

	return;
}


/*------------------------------------------------------------------------------
GEOMETRIC_K - space frame geometric stiffness matrix, global coordnates 20dec07
------------------------------------------------------------------------------*/
void geometric_K(
	float **k, float *x, float *y, float *z, float *r
	, float L, float Le
	, int j1, int j2, float Ax, float Asy, float Asz, float J
	, float Iy, float Iz, float E, float G, float p, float T
	, int shear
){
	float t1, t2, t3, t4, t5, t6, t7, t8, t9; /* coord Xformn */
	float **kg;
	float Ksy,Ksz,Dsy,Dsz; /* shear deformation coefficients	*/
	int i, j;

	coord_trans(
		x, y, z, L, j1, j2, &t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p 
	);

	kg = matrix(1,12,1,12);
	for (i=1;i<=12;i++)	for (j=1;j<=12;j++)	kg[i][j] = 0.0;

	if ( shear ) {
		Ksy = 12.*E*Iz / (G*Asy*Le*Le);
		Ksz = 12.*E*Iy / (G*Asz*Le*Le);
		Dsy = (1+Ksy)*(1+Ksy);
		Dsz = (1+Ksz)*(1+Ksz);
	} else{
		Ksy = Ksz = 0.0;
		Dsy = Dsz = 1.0;
	}


        kg[1][1]  = kg[7][7]   = T/L;

	kg[2][2]  = kg[8][8]   = T/L*(1.2+2.0*Ksy+Ksy*Ksy)/Dsy;
	kg[3][3]  = kg[9][9]   = T/L*(1.2+2.0*Ksz+Ksz*Ksz)/Dsz;
	kg[4][4]  = kg[10][10] = T/L*J/Ax;
	kg[5][5]  = kg[11][11] = T*L*(2.0/15.0+Ksz/6.0+Ksz*Ksz/12.0)/Dsz;
	kg[6][6]  = kg[12][12] = T*L*(2.0/15.0+Ksy/6.0+Ksy*Ksy/12.0)/Dsy;

        kg[1][7]  = kg[7][1]   = -T/L;
 
	kg[5][3]  = kg[3][5]   =  kg[11][3] = kg[3][11] = -T/10.0/Dsz;
	kg[9][5]  = kg[5][9]   =  kg[11][9] = kg[9][11] =  T/10.0/Dsz;
	kg[6][2]  = kg[2][6]   =  kg[12][2] = kg[2][12] =  T/10.0/Dsy;
	kg[8][6]  = kg[6][8]   =  kg[12][8] = kg[8][12] = -T/10.0/Dsy;

        kg[4][10] = kg[10][4]  = -kg[4][4];

	kg[8][2]  = kg[2][8]   = -T/L*(1.2+2.0*Ksy+Ksy*Ksy)/Dsy;
	kg[9][3]  = kg[3][9]   = -T/L*(1.2+2.0*Ksz+Ksz*Ksz)/Dsz;

	kg[11][5] = kg[5][11]  = -T*L*(1.0/30.0+Ksz/6.0+Ksz*Ksz/12.0)/Dsz;
	kg[12][6] = kg[6][12]  = -T*L*(1.0/30.0+Ksy/6.0+Ksy*Ksy/12.0)/Dsy;

/*	save_matrix ( 12, 12, kg, "kg" ); /* element geom. stiffness matrix */

	atma ( t1,t2,t3,t4,t5,t6,t7,t8,t9, kg, r[j1],r[j2] );	/* globalize */

	/* check and enforce symmetry */ 

	for (i=1; i<=12; i++)
	    for (j=i+1; j<=12; j++) 
		if ( kg[i][j] != kg[j][i] ) {
			kg[i][j] = kg[j][i] = 0.5 * ( kg[i][j] + kg[j][i] );
/*			fprintf(stderr,"kg[%d][%d] = %e    ",i,j,kg[i][j] ); */
/*			fprintf(stderr,"kg[%d][%d] = %e  \n",j,i,kg[j][i] ); */
		}

/*	save_matrix ( 12, 12, kg, "kgt" );   /* transformed element matx */

	/* add geometric stiffness matrix to elastic stiffness matrix ... */

	for (i=1; i<=12; i++)   for (j=1; j<=12; j++)	k[i][j] += kg[i][j];

	free_matrix(kg,1,12,1,12);

	return;
}


/*------------------------------------------------------------------------------
END_RELEASE - apply matrix condensation for one member end force release 20nov04
------------------------------------------------------------------------------*/
void end_release ( X, r )
float	**X;
int	r;
{
	int	i,j;

	for (i=1; i<=12; i++)  if ( i != r) 
		for (j=1; j<=12; j++) if ( j != r) 
			X[i][j]  =  X[i][j]  -  X[i][r] * X[r][j] / X[r][r];
	
	for (i=1; i<=12; i++) 
		X[r][i] = X[i][r] = 0.0;

	return;
}


/*------------------------------------------------------------------------------
APPLY_REACTIONS -  apply fixed joint displacement boundary conditions	23feb94
The original external load vector, Fo is returned unchanged;
The load vector modified for prescribed displacements, Dp, is returned as F
------------------------------------------------------------------------------*/
void apply_reactions ( DoF, R, Dp, Fo, F, K )
int	DoF, *R;
float	*Dp, *Fo, *F, **K;
{
	int	i,j;

	for (i=1; i<=DoF; i++)  F[i] = Fo[i];

	for (i=1; i<=DoF; i++) {		/* modify the force vector */
		if ( R[i] ) {
			F[i] = Dp[i];
		} else {
			for (j=1; j<=DoF; j++) 
				if ( R[j] ) F[i] -= K[i][j]*Dp[j];
		}
	}
			
	for (i=1; i<=DoF; i++) {	/*  modify the stiffness matrix  */
		if ( R[i] ) {
			for (j=1; j<=DoF; j++) {
				if ( i == j )	K[i][j] = 1.0;
				else		K[i][j] = K[j][i] = 0.0;
			}
		}
	}
	return;
}


/*----------------------------------------------------------------------------
SOLVE_SYSTEM  -  solve {F} =   [K]{D} via L D L' decomposition        27dec01
----------------------------------------------------------------------------*/
void solve_system( K, D, F, DoF, ok )
float	**K, *D, *F; 
int	DoF, *ok;
{
	float	*diag,		/* diagonal vector of the L D L' decomp. */
		error=1.0;	/* error in the solution		*/

	diag = vector ( 1, DoF );

	ldl_dcmp( K, DoF, diag, F, D, 1, 0, ok );	/*  L D L'  decomp */
	if ( *ok < 0 ) {
	 	fprintf(stderr," Make sure that all six");
		fprintf(stderr," rigid body translations are restrained!\n");
		/* exit(1); */
	} else {				/* back substitute for D */
		ldl_dcmp( K, DoF, diag, F, D, 0, 1, ok ); /* LDL'  back-sub */
	        fprintf(stderr,"    LDL' RMS matrix precision:");
		error = *ok = 1;
		do {					/* improve solution */
			ldl_mprove ( K, DoF, diag, F, D, &error, ok );
			fprintf(stderr,"%9.2e", error );
		} while ( *ok );
	        fprintf(stderr,"\n");
	}

	free_vector( diag, 1, DoF );

	return;
}


/*------------------------------------------------------------------------------
END_FORCES  -  evaluate the member end forces for every member		23feb94
------------------------------------------------------------------------------*/
void end_forces ( Q, nM, x,y,z, L, Le,
			J1,J2, Ax, Asy,Asz, J, Iy,Iz, E, G, p, D, shear, geom )
int	nM, *J1, *J2, shear, geom;
float	**Q, *x, *y, *z, *L, *Le, *Ax, *Asy, *Asz, *J, *Iy, *Iz, *E, *G, *p, *D;
{
	float	*s;
	int	i,j;

	s = vector(1,12);

	for (i=1; i <= nM; i++) {

     		member_force ( s, i, x,y,z, L[i], Le[i], J1[i], J2[i],
				Ax[i], Asy[i], Asz[i], J[i], Iy[i], Iz[i],
					E[i], G[i], p[i], D, shear, geom );

		for (j=1; j<=12; j++)	Q[i][j] = s[j];
	}
	free_vector(s,1,12);
	return;
}


/*------------------------------------------------------------------------------
MEMBER_FORCE  -  evaluate the end forces for a member			12nov02
------------------------------------------------------------------------------*/
void member_force(
	float *s, int M, float *x, float *y, float *z, float L, float Le
	, int j1, int j2, float Ax, float Asy, float Asz, float J
	, float Iy, float Iz, float E, float G, float p, float *D
	, int shear, int geom
){
	float	t1, t2, t3, t4, t5, t6, t7, t8, t9, /* coord Xformn	*/
		d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12,
		x1, y1, z1, x2, y2, z2,	/* joint coordinates	*/
		Ls,			/* stretched length of element */
		Ksy, Ksz, Dsy, Dsz,	/* shear deformation coeff's	*/
		T;		/* normal force for geometric stiffness */

	coord_trans ( x, y, z, L, j1, j2,
			&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p );

	x1 = x[j1];	y1 = y[j1];	z1 = z[j1];
	x2 = x[j2];	y2 = y[j2];	z2 = z[j2];

	j1 = 6*(j1-1);	j2 = 6*(j2-1);

	d1  = D[j1+1];	d2  = D[j1+2];	d3  = D[j1+3];
	d4  = D[j1+4];	d5  = D[j1+5];	d6  = D[j1+6];
	d7  = D[j2+1];	d8  = D[j2+2];	d9  = D[j2+3];
	d10 = D[j2+4];	d11 = D[j2+5];	d12 = D[j2+6];

	if ( shear ) {
		Ksy = 12.*E*Iz / (G*Asy*Le*Le);
		Ksz = 12.*E*Iy / (G*Asz*Le*Le);
		Dsy = (1+Ksy)*(1+Ksy);
		Dsz = (1+Ksz)*(1+Ksz);
	} else {
		Ksy = Ksz = 0.0;
		Dsy = Dsz = 1.0;
	}

	if ( geom ) {
		Ls    = pow((x2+d7-x1-d1),2.0) + 
			pow((y2+d8-y1-d2),2.0) + 
			pow((z2+d9-z1-d3),2.0);
		Ls = sqrt(Ls) + Le - L;
		T  = Ax*E*log(Ls/Le);	/* true strain	*/
		T  =  (Ax*E/Le)*( (d7-d1)*t1 + (d8-d2)*t2 + (d9-d3)*t3 ); 
	}
	else        T =  0.0;

	if ( geom )
		s[1] = -T;
	else
		s[1]  =  -(Ax*E/Le)*( (d7-d1)*t1 + (d8-d2)*t2 + (d9-d3)*t3 );
	s[2]  = -( 12.*E*Iz/(Le*Le*Le*(1.+Ksy)) + 
		   T/L*(1.2+2.0*Ksy+Ksy*Ksy)/Dsy ) *
				( (d7-d1)*t4 + (d8-d2)*t5 + (d9-d3)*t6 )
		+ (6.*E*Iz/(Le*Le*(1.+Ksy)) + T/10.0/Dsy) *
				( (d4+d10)*t7 + (d5+d11)*t8 + (d6+d12)*t9 );
	s[3]  = -(12.*E*Iy/(Le*Le*Le*(1.+Ksz)) + 
		  T/L*(1.2+2.0*Ksz+Ksz*Ksz)/Dsz ) * 
				( (d7-d1)*t7  + (d8-d2)*t8  + (d9-d3)*t9 )
		- ( 6.*E*Iy/(Le*Le*(1.+Ksz)) + T/10.0/Dsz ) *
				( (d4+d10)*t4 + (d5+d11)*t5 + (d6+d12)*t6 );
	s[4]  =   -(G*J/Le) * ( (d10-d4)*t1 + (d11-d5)*t2 + (d12-d6)*t3 );
	s[5]  =   (6.*E*Iy/(Le*Le*(1.+Ksz)) + T/10.0/Dsz ) * 
				( (d7-d1)*t7 + (d8-d2)*t8 + (d9-d3)*t9 )
		+ ( (4.+Ksz)*E*Iy/(Le*(1.+Ksz)) + 
		    T*L*(2.0/15.0+Ksz/6.0+Ksz*Ksz/12.0)/Dsz ) *
				(d4 *t4 + d5 *t5 + d6 *t6)
		+ ((2.-Ksz)*E*Iy/(Le*(1.+Ksz)) -  
		    T*L*(1.0/30.0+Ksz/6.0+Ksz*Ksz/12.0)/Dsz ) *
				(d10*t4 + d11*t5 + d12*t6);
	s[6]  =  -( 6.*E*Iz/(Le*Le*(1.+Ksy)) + T/10.0/Dsy ) *
				( (d7-d1)*t4 + (d8-d2)*t5 + (d9-d3)*t6 )
		+ ((4.+Ksy)*E*Iz/(Le*(1.+Ksy)) + 
		    T*L*(2.0/15.0+Ksy/6.0+Ksy*Ksy/12.0)/Dsy ) *
				( d4 *t7 + d5 *t8 + d6 *t9 )
		+ ((2.-Ksy)*E*Iz/(Le*(1.+Ksy)) - 
		    T*L*(1.0/30.0+Ksy/6.0+Ksy*Ksy/12.0)/Dsy ) * 
				( d10*t7 + d11*t8 + d12*t9 );
	s[7]  = -s[1];
	s[8]  = -s[2]; 
	s[9]  = -s[3]; 
	s[10] = -s[4]; 

	s[11] =   ( 6.*E*Iy/(Le*Le*(1.+Ksz)) + T/10.0/Dsz ) * 
				( (d7-d1)*t7 + (d8-d2)*t8 + (d9-d3)*t9 )
		+ ((4.+Ksz)*E*Iy/(Le*(1.+Ksz)) + 
		    T*L*(2.0/15.0+Ksz/6.0+Ksz*Ksz/12.0)/Dsz ) * 
				( d10*t4 + d11*t5 + d12*t6 )
		+ ((2.-Ksz)*E*Iy/(Le*(1.+Ksz)) - 
		    T*L*(1.0/30.0+Ksz/6.0+Ksz*Ksz/12.0)/Dsz ) *
				( d4 *t4 + d5 *t5 + d6 *t6 );
	s[12] =  -(6.*E*Iz/(Le*Le*(1.+Ksy)) + T/10.0/Dsy ) *
				( (d7-d1)*t4 + (d8-d2)*t5 + (d9-d3)*t6 )
		+ ((4.+Ksy)*E*Iz/(Le*(1.+Ksy)) + 
		    T*L*(2.0/15.0+Ksy/6.0+Ksy*Ksy/12.0)/Dsy ) * 
				( d10*t7 + d11*t8 + d12*t9 )
		+ ((2.-Ksy)*E*Iz/(Le*(1.+Ksy)) - 
		    T*L*(1.0/30.0+Ksy/6.0+Ksy*Ksy/12.0)/Dsy ) * 
				( d4 *t7 + d5 *t8 + d6 *t9 );

	return;
}


/*----------------------------------------------------------------------------- 
EQUILIBRIUM  -  perform an equilibrium check, F returned as reactions   18sep02
------------------------------------------------------------------------------*/
void equilibrium ( x, y, z, L, J1, J2, F, R, p, Q, feF, nM, DoF, err )
int	*J1, *J2, *R, nM, DoF;
float	*x, *y, *z, *L, *F, **Q, *p, **feF, *err;
{
	float   t1, t2, t3, t4, t5, t6, t7, t8, t9,	/* 3D coord Xformn */
		f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12,
		den = 0.0;			/* normalizes the error	*/
	int	i,j,m, j1, j2;

	for (j=1; j<=DoF; j++)	if ( R[j] == 0 )   den += ( F[j]*F[j] );
	den = sqrt ( den / (float) DoF );
	if ( den <= 0 ) den = 1;

	for (m=1; m <= nM; m++) {	/* loop over all members */

		j1 = J1[m];	j2 = J2[m];

		coord_trans ( x, y, z, L[m], j1, j2,
			&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p[m] );

		j1 = 6*(j1-1);	j2 = 6*(j2-1);

	/* subtract internal forces from external loads	*/

							/* {F} = [T]' {Q} */
		F[j1+1] -= ( Q[m][1] *t1 + Q[m][2] *t4 + Q[m][3] *t7 );
		F[j1+2] -= ( Q[m][1] *t2 + Q[m][2] *t5 + Q[m][3] *t8 );
		F[j1+3] -= ( Q[m][1] *t3 + Q[m][2] *t6 + Q[m][3] *t9 );
		F[j1+4] -= ( Q[m][4] *t1 + Q[m][5] *t4 + Q[m][6] *t7 );
		F[j1+5] -= ( Q[m][4] *t2 + Q[m][5] *t5 + Q[m][6] *t8 );
		F[j1+6] -= ( Q[m][4] *t3 + Q[m][5] *t6 + Q[m][6] *t9 );

		F[j2+1] -= ( Q[m][7] *t1 + Q[m][8] *t4 + Q[m][9] *t7 );
		F[j2+2] -= ( Q[m][7] *t2 + Q[m][8] *t5 + Q[m][9] *t8 );
		F[j2+3] -= ( Q[m][7] *t3 + Q[m][8] *t6 + Q[m][9] *t9 );
		F[j2+4] -= ( Q[m][10]*t1 + Q[m][11]*t4 + Q[m][12]*t7 );
		F[j2+5] -= ( Q[m][10]*t2 + Q[m][11]*t5 + Q[m][12]*t8 );
		F[j2+6] -= ( Q[m][10]*t3 + Q[m][11]*t6 + Q[m][12]*t9 );

	/* subtract fixed end forces (equivalent loads) from internal loads */
	
		f1 = feF[m][1];		f2 = feF[m][2];		f3 = feF[m][3];
		f4 = feF[m][4];		f5 = feF[m][5];		f6 = feF[m][6];
		f7 = feF[m][7];		f8 = feF[m][8];		f9 = feF[m][9];
		f10= feF[m][10];	f11= feF[m][11];	f12= feF[m][12];

		Q[m][1]  -= ( f1 *t1 + f2 *t2 + f3 *t3 );     /* {Q} = [T]{F} */
		Q[m][2]  -= ( f1 *t4 + f2 *t5 + f3 *t6 );
		Q[m][3]  -= ( f1 *t7 + f2 *t8 + f3 *t9 );
		Q[m][4]  -= ( f4 *t1 + f5 *t2 + f6 *t3 );
		Q[m][5]  -= ( f4 *t4 + f5 *t5 + f6 *t6 );
		Q[m][6]  -= ( f4 *t7 + f5 *t8 + f6 *t9 );

		Q[m][7]  -= ( f7 *t1 + f8 *t2 + f9 *t3 );
		Q[m][8]  -= ( f7 *t4 + f8 *t5 + f9 *t6 );
		Q[m][9]  -= ( f7 *t7 + f8 *t8 + f9 *t9 );
		Q[m][10] -= ( f10*t1 + f11*t2 + f12*t3 );
		Q[m][11] -= ( f10*t4 + f11*t5 + f12*t6 );
		Q[m][12] -= ( f10*t7 + f11*t8 + f12*t9 );


	}

	*err = 0.0;
	for (j=1; j<=DoF; j++)	if ( R[j] == 0 )	*err += ( F[j]*F[j] );
	*err = sqrt ( *err / (float) DoF ) / den;
	fprintf(stderr,"  RMS relative equilibrium precision: %9.3e\n", *err );

	return;
}

/*------------------------------------------------------------------------------
ASSEMBLE_M  -  assemble global mass matrix from element mass & inertia  24nov98
------------------------------------------------------------------------------*/
void assemble_M ( M, DoF, nJ, nM, x,y,z,r, L, J1,J2, Ax, J,Iy,Iz, p,
					d, BMs, JMs, JMx, JMy, JMz, lump )
int     DoF, nM;
float   **M, *x,*y,*z,*r, *L, *Ax, *J, *Iy,*Iz, *p, *d, *BMs, *JMs, *JMx,*JMy,*JMz;
int     *J1, *J2, lump;
{
	float   **mass,	    /* element mass matrix in global coord */
		**matrix();
	int     **ind,	  /* member-structure DoF index table     */
		**imatrix(),
		i, j, ii, jj, l, ll, m;

	for (i=1; i<=DoF; i++)  for (j=1; j<=DoF; j++)  M[i][j] = 0.0;

	mass   =  matrix(1,12,1,12);
	ind    = imatrix(1,12,1,nM);


	for ( i=1; i<= nM; i++ ) {
		ind[1][i] = 6*J1[i] - 5;	ind[7][i]  = 6*J2[i] - 5;
		ind[2][i] = ind[1][i] + 1;      ind[8][i]  = ind[7][i] + 1;
		ind[3][i] = ind[1][i] + 2;      ind[9][i]  = ind[7][i] + 2;
		ind[4][i] = ind[1][i] + 3;      ind[10][i] = ind[7][i] + 3;
		ind[5][i] = ind[1][i] + 4;      ind[11][i] = ind[7][i] + 4;
		ind[6][i] = ind[1][i] + 5;      ind[12][i] = ind[7][i] + 5;
	}

	for ( m = 1; m <= nM; m++ ) {

		if ( lump )	lumped_M ( mass, x,y,z, L[m], J1[m], J2[m],
				Ax[m], J[m], Iy[m], Iz[m], d[m], BMs[m], p[m]);
		else		consistent_M ( mass, x,y,z,r,L[m], J1[m], J2[m],
				Ax[m], J[m], Iy[m], Iz[m], d[m], BMs[m], p[m]);

		for ( l=1; l <= 12; l++ ) {
			ii = ind[l][m];
			for ( ll=1; ll <= 12; ll++ ) {
				jj = ind[ll][m];
				M[ii][jj] += mass[l][ll];
			}
		}
	}
	for ( j = 1; j <= nJ; j++ ) {
		i = 6*(j-1);
		M[i+1][i+1] += JMs[j];
		M[i+2][i+2] += JMs[j];
		M[i+3][i+3] += JMs[j];
		M[i+4][i+4] += JMx[j];
		M[i+5][i+5] += JMy[j];
		M[i+6][i+6] += JMz[j];
	}

	for (i=1; i<= DoF; i++) {
		if ( M[i][i] <= 0.0 ) {
			fprintf(stderr,"  error: Non pos-def mass matrix\n");
			fprintf(stderr,"  M[%d][%d] = %lf\n", i,i, M[i][i] );
		}
	}
	free_matrix ( mass,1,12,1,12);
	free_imatrix( ind,1,12,1,nM);
	return;
}


/*------------------------------------------------------------------------------
LUMPED_M  -  space frame element lumped mass matrix in global coordnates 7apr94
------------------------------------------------------------------------------*/
void lumped_M(
	float **m, float *x, float *y, float *z
	, float L, int j1, int j2, float Ax, float J, float Iy, float Iz
	, float d, float p, float BMs
){
	float   t1, t2, t3, t4, t5, t6, t7, t8, t9,     /* coord Xformn */
		t, ry,rz, po;	/* translational, rotational & polar inertia */
	int     i, j;

	coord_trans ( x, y, z, L, j1, j2,
				&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p );

			/* rotatory inertia of extra mass is neglected */

	t = ( d*Ax*L + BMs ) / 2.0;
	ry = d*Iy*L / 2.0;
	rz = d*Iz*L / 2.0;
	po = d*L*J / 2.0;		/* assumes simple cross-section	*/

	for (i=1;i<=12;i++)	for (j=1;j<=12;j++)	m[i][j] = 0.0;

	m[1][1] = m[2][2] = m[3][3] = m[7][7] = m[8][8] = m[9][9] = t;

	m[4][4] = m[10][10] = po*t1*t1 + ry*t4*t4 + rz*t7*t7;
	m[5][5] = m[11][11] = po*t2*t2 + ry*t5*t5 + rz*t8*t8;
	m[6][6] = m[12][12] = po*t3*t3 + ry*t6*t6 + rz*t9*t9;

	m[4][5] = m[5][4] = m[10][11] = m[11][10] =po*t1*t2 +ry*t4*t5 +rz*t7*t8;
	m[4][6] = m[6][4] = m[10][12] = m[12][10] =po*t1*t3 +ry*t4*t6 +rz*t7*t9;
	m[5][6] = m[6][5] = m[11][12] = m[12][11] =po*t2*t3 +ry*t5*t6 +rz*t8*t9;

	return;
}


/*------------------------------------------------------------------------------
CONSISTENT_M  -  space frame consistent mass matrix in global coordnates 2oct97
		 does not include shear deformations
------------------------------------------------------------------------------*/
void consistent_M(
	float **m, float *x, float *y, float *z, float *r, float L
	, int j1, int j2, float Ax, float J, float Iy, float Iz, float d
	, float BMs, float p
){
	float   t1, t2, t3, t4, t5, t6, t7, t8, t9,     /* coord Xformn */
		t, ry, rz, po;	/* translational, rotational & polar inertia */
	int     i, j;

	coord_trans ( x, y, z, L, j1, j2,
				&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p );

	t  =  d*Ax*L;	
	ry =  d*Iy;
	rz =  d*Iz;
	po =  d*J*L;

	for (i=1;i<=12;i++)	for (j=1;j<=12;j++)	m[i][j] = 0.0;

	m[1][1]  = m[7][7]   = t/3.;
	m[2][2]  = m[8][8]   = 13.*t/35. + 6.*rz/(5.*L);
	m[3][3]  = m[9][9]   = 13.*t/35. + 6.*ry/(5.*L);
	m[4][4]  = m[10][10] = po/3.;
	m[5][5]  = m[11][11] = t*L*L/105. + 2.*L*ry/15.;
	m[6][6]  = m[12][12] = t*L*L/105. + 2.*L*rz/15.;

	m[5][3]  = m[3][5]   = -11.*t*L/210. - ry/10.;
	m[6][2]  = m[2][6]   =  11.*t*L/210. + rz/10.;
	m[7][1]  = m[1][7]   =  t/6.;

	m[8][6]  = m[6][8]   =  13.*t*L/420. - rz/10.;
	m[9][5]  = m[5][9]   = -13.*t*L/420. + ry/10.;
	m[10][4] = m[4][10]  =  po/6.; 
	m[11][3] = m[3][11]  =  13.*t*L/420. - ry/10.;
	m[12][2] = m[2][12]  = -13.*t*L/420. + rz/10.;

	m[11][9] = m[9][11]  =  11.*t*L/210. + ry/10.;
	m[12][8] = m[8][12]  = -11.*t*L/210. - rz/10.;

	m[8][2]  = m[2][8]   =  9.*t/70. - 6.*rz/(5.*L);
	m[9][3]  = m[3][9]   =  9.*t/70. - 6.*ry/(5.*L);
	m[11][5] = m[5][11]  = -L*L*t/140. - ry*L/30.;
	m[12][6] = m[6][12]  = -L*L*t/140. - rz*L/30.;

/*	save_matrix ( 12, 12, m, "mo" );	/* element mass matrix */

	atma ( t1,t2,t3,t4,t5,t6,t7,t8,t9, m, r[j1],r[j2] );	/* globalize */

	/* check and enforce symmetry */ 

	for (i=1; i<=12; i++)
	    for (j=i+1; j<=12; j++)
		if ( m[i][j] != m[j][i] ) {
			m[i][j] = m[j][i] = 0.5 * ( m[i][j] + m[j][i] );
/*			fprintf(stderr,"m[%d][%d] = %e    ",i,j,m[i][j] ); */
/*			fprintf(stderr,"m[%d][%d] = %e  \n",j,i,m[j][i] ); */
		}

			/* rotatory inertia of extra mass is neglected */

	for (i=1; i<=3; i++)	m[i][i] += BMs/2.;
	for (i=7; i<=9; i++)	m[i][i] += BMs/2.;

/*	save_matrix ( 12, 12, m, "mt" );	/* transformed matrix */

	return;
}

/*------------------------------------------------------------------------------
ATMA  -  perform the coordinate transformation from local to global 	6jan96
         include effects of a finite joint radii, r1 and r2.            9dec04
------------------------------------------------------------------------------*/
void atma ( t1, t2, t3, t4, t5, t6, t7, t8, t9, m, r1, r2 )
float	t1, t2, t3, t4, t5, t6, t7, t8, t9, **m, r1, r2;
{
	float	**a, **ma, **matrix();
	int	i,j,k;

	a  = matrix(1,12,1,12);
	ma = matrix(1,12,1,12);

	for (i=1; i<=12; i++)
	    for (j=i; j<=12; j++) 
		ma[j][i] = ma[i][j] = a[j][i] = a[i][j] = 0.0;

	for (i=0; i<=3; i++) {
		a[3*i+1][3*i+1] = t1;
		a[3*i+1][3*i+2] = t2;
		a[3*i+1][3*i+3] = t3;
		a[3*i+2][3*i+1] = t4;
		a[3*i+2][3*i+2] = t5;
		a[3*i+2][3*i+3] = t6;
		a[3*i+3][3*i+1] = t7;
		a[3*i+3][3*i+2] = t8;
		a[3*i+3][3*i+3] = t9;
	}

/*
	a[5][1] =  r1*t7; 
	a[5][2] =  r1*t8; 
	a[5][3] =  r1*t9; 
	a[6][1] = -r1*t4; 
	a[6][2] = -r1*t5; 
	a[6][3] = -r1*t6; 

	a[11][7] = -r2*t7; 
	a[11][8] = -r2*t8; 
	a[11][9] = -r2*t9; 
	a[12][7] =  r2*t4; 
	a[12][8] =  r2*t5; 
	a[12][9] =  r2*t6; 
*/


/*	save_matrix( 12, 12, a, "aa");		/* save cord xfmtn */

	for (j=1; j <= 12; j++)				/* MT = M T	*/
	    for (i=1; i <= 12; i++)
		for (k=1; k <= 12; k++)    ma[i][j] += m[i][k] * a[k][j];

/*	save_matrix( 12, 12, ma, "ma");	/* partial transformation */

	for (i=1; i<=12; i++)   for (j=i; j<=12; j++)   m[j][i] = m[i][j] = 0.0;

	for (j=1; j <= 12; j++)				/* T'MT = T' MT	*/
	    for (i=1; i <= 12; i++)
		for (k=1; k <= 12; k++)    m[i][j] += a[k][i] * ma[k][j];

/*	save_matrix( 12, 12, m, "atma");		/* debug atma */

	free_matrix(a, 1,12,1,12);
	free_matrix(ma,1,12,1,12);
	return;
} 


/*------------------------------------------------------------------------------
CONDENSE - static condensation of stiffness matrix from NxN to nxn    30aug01
------------------------------------------------------------------------------*/
void condense ( A, N, q, n, Ac)
float	**A, **Ac;
int	N, n, *q;
{
	float	**Arr, **Arq;
	int	i,j,k, ri,rj,qi,qj, ok, 
		*r;

	r    = ivector(1,N-n);
	Arr  = matrix(1,N-n,1,N-n);
	Arq  = matrix(1,N-n,1,n);

	k = 1;
	for (i=1; i<=N; i++) {
		ok = 1;
		for (j=1; j<=n; j++) {
			if ( q[j] == i ) {
				ok = 0;
				break;
			}
		}
		if ( ok )	r[k++] = i;
	}

	for (i=1; i<=N-n; i++) {
		for (j=i; j<=N-n; j++) { /* use only upper triangle of A */
			ri = r[i];
			rj = r[j];
			if ( ri <= rj )	Arr[j][i] = Arr[i][j] = A[ri][rj];
		}
	}

	for (i=1; i<=N-n; i++) {
		for (j=1; j<=n; j++) {	/* use only upper triangle of A */
			ri = r[i];
			qj = q[j];
			if ( ri < qj )	Arq[i][j] = A[ri][qj];
			else		Arq[i][j] = A[qj][ri];
		}
	}

	xtaiy( Arq, Arr, Arq, N-n, n, Ac );

	for (i=1; i<=n; i++) {
		for (j=i; j<=n; j++) { /* use only upper triangle of A */
			qi = q[i];
			qj = q[j];
			if ( qi <= qj ) Ac[j][i]=Ac[i][j] = A[qi][qj]-Ac[i][j];
		}
	}

	free_ivector ( r,   1,N-n );
	free_matrix  ( Arr, 1,N-n,1,N-n );
	free_matrix  ( Arq, 1,N-n,1,n );

	return;
}


/*---------------------------------------------------------------------------- 
GUYAN  -   generalized Guyan reduction of mass and stiffness matrices    6jun07
           matches the response at a particular frequency, sqrt(L)/2/pi
           Guyan, Robert J., ``Reduction of Stiffness and Mass Matrices,''
           AIAA Journal, Vol. 3, No. 2 (1965) p 380.
-----------------------------------------------------------------------------*/
void guyan ( M, K, N, q, n, Mc, Kc, w2 )
float	**M, **K, **Mc, **Kc, w2;
int	N, n, *q;
{
	float	**Drr, **Drq, **invDrrDrq, **T;
	int	i,j,k, ri,rj,qj, ok, 
		*r;

	r   = ivector(1,N-n);
	Drr =  matrix(1,N-n,1,N-n);
	Drq =  matrix(1,N-n,1,n);
	invDrrDrq = matrix(1,N-n,1,n);	/* inv(Drr) * Drq	*/
	T   = matrix(1,N,1,n);	/* coordinate transformation matrix	*/

	w2 = 4.0 * PI * PI * w2 * w2;	/* eigen-value ... omega^2 	*/

	/* find "remaining" (r) degrees of freedom, not "qondensed" (q)	*/
	k = 1;
	for (i=1; i<=N; i++) {
		ok = 1;
		for (j=1; j<=n; j++) {
			if ( q[j] == i ) {
				ok = 0;
				break;
			}
		}
		if ( ok )	r[k++] = i;
	}

	for (i=1; i<=N-n; i++) {
		for (j=1; j<=N-n; j++) { /* use only upper triangle of K,M */
			ri = r[i];
			rj = r[j];
			if ( ri <= rj )	
				Drr[j][i] = Drr[i][j] = K[ri][rj]-w2*M[ri][rj];
			else	Drr[j][i] = Drr[i][j] = K[rj][ri]-w2*M[rj][ri];
		}
	}

	for (i=1; i<=N-n; i++) {
		for (j=1; j<=n; j++) {	/* use only upper triangle of K,M */
			ri = r[i];
			qj = q[j];
			if ( ri < qj )	Drq[i][j] = K[ri][qj] - w2*M[ri][qj];
			else		Drq[i][j] = K[qj][ri] - w2*M[qj][ri];
		}
	}

	invAB(Drr, Drq, N-n, n, invDrrDrq, &ok);	/* inv(Drr) * Drq	*/

	/* coordinate transformation matrix	*/	
	for (i=1; i<=n; i++) {
		for (j=1; j<=n; j++)	T[q[i]][j] =  0.0;
		T[q[i]][i] = 1.0;
	}	
	for (i=1; i<=N-n; i++) 
		for (j=1; j<=n; j++)	T[r[i]][j] = -invDrrDrq[i][j];

	xtAx ( K, T, Kc, N, n );		/* Kc = T' * K * T	*/

	xtAx ( M, T, Mc, N, n );		/* Mc = T' * M * T	*/

	free_ivector ( r,   1, N-n );
	free_matrix  ( Drr, 1,N-n,1,N-n );
	free_matrix  ( Drq, 1,N-n,1,n );
	free_matrix  ( invDrrDrq, 1,N-n,1,N-n );
	free_matrix  ( T, 1,N-n,1,n );

	return;
}


/*---------------------------------------------------------------------------- 
DYN_CONDEN - dynamic condensation of mass and stiffness matrices    8oct01
	     matches the response at a set of frequencies
WARNING: Kc and Mc may be ill-conditioned, and possibly non-positive def.
-----------------------------------------------------------------------------*/
void dyn_conden ( M, K, N, R, p, n, Mc, Kc, V, f, m )
float	**M, **K, **Mc, **Kc, **V, *f;
int	*R, N, n, *p, *m;
{
	float	**P, **invP,
		traceM = 0, traceMc = 0, 
		Aij;		/* temporary storage for matrix mult. */
	int	i,j,k, pi, ok; 

	P    =  matrix(1,n,1,n);
	invP =  matrix(1,n,1,n);

	for (i=1; i<=n; i++)	/* first n modal vectors at primary DoF's */
		for (j=1; j<=n; j++)
			P[i][j] = V[p[i]][m[j]];

	pseudo_inv ( P, invP, n, n, 1e-9 );

	for (i=1; i<=N; i++) if ( !R[i] ) traceM += M[i][i];

	for (i=1; i<=n; i++) { 		/* compute inv(P)' * I * inv(P)	*/
	    for (j=1; j<=n; j++) {
		Aij = 0.0;
	        for (k=1; k<=n; k++)
			Aij += invP[k][i] * invP[k][j];
		Mc[i][j] = Aij;
	    }
	}

	for (i=1; i<=n; i++) traceMc += Mc[i][i];

	for (i=1; i<=n; i++) { 		/* compute inv(P)' * W^2 * inv(P) */
	    for (j=1; j<=n; j++) {
		Aij = 0.0;
	        for (k=1; k<=n; k++) 
		    Aij += invP[k][i] * 4.0*PI*PI*f[m[k]]*f[m[k]] * invP[k][j];
		Kc[i][j] = Aij;
	    }
	}

	for (i=1; i<=n; i++)
	       for (j=1; j<=n; j++)
		       Mc[i][j] *= (traceM / traceMc);

	for (i=1; i<=n; i++)
	       for (j=1; j<=n; j++)
		       Kc[i][j] *= (traceM / traceMc);

	free_matrix  ( P,    1,n,1,n);
	free_matrix  ( invP, 1,n,1,n);

	return;
}


/*------------------------------------------------------------------------------
INVAB  -  calculate product inv(A) * B  
        A is n by n      B is n by m                                    6jun07
------------------------------------------------------------------------------*/
void invAB(
	float **A, float **B
	, int n, int m, float **AiB
	, int *ok
){
	float	*diag, *b, *x, error;
	int	i,j,k, disp=1;

	diag = vector(1,n);
	x    = vector(1,n);
	b    = vector(1,n);

	for (i=1; i<=n; i++) diag[i] = x[i] = 0.0;

	ldl_dcmp( A, n, diag, b, x, 1, 0, ok );		/*  L D L'  decomp */
	if ( *ok < 0 ) {
	 	fprintf(stderr," Make sure that all six");
		fprintf(stderr," rigid body translations are restrained!\n");
	}

	for (j=1; j<=m; j++) {

		for (k=1; k<=n; k++)  b[k] = B[k][j];
		ldl_dcmp( A, n, diag, b, x, 0, 1, ok ); /*  L D L'  bksbtn */

		if (disp)
		 fprintf(stderr,"    LDL' RMS matrix precision:");
		error = *ok = 1;
		do {					/*improve the solution*/
			ldl_mprove ( A, n, diag, b, x, &error, ok );
			if (disp) fprintf(stderr,"%9.2e", error );
		} while ( *ok );
		if (disp) fprintf(stderr,"\n");

		for (i=1; i<=n; i++)	AiB[i][j] = x[i];
	}

	free_vector(diag,1,n);
	free_vector(x,1,n);
	free_vector(b,1,n);
	return;
}


/*------------------------------------------------------------------------------
XTAIY  -  calculate quadratic form with inverse matrix   X' * inv(A) * Y  
          A is n by n    X is n by m     Y is n by m                    15sep01
------------------------------------------------------------------------------*/
void xtaiy(
	float **X, float **A, float **Y
	, int n, int m, float **Ac
){
	float	*diag, *x, *y, error;
	int	i,j,k, ok, disp=0;

	diag = vector(1,n);
	x    = vector(1,n);
	y    = vector(1,n);

	for (i=1; i<=n; i++) diag[i] = x[i] = 0.0;

	ldl_dcmp( A, n, diag, y, x, 1, 0, &ok );	/*  L D L'  decomp */

	for (j=1; j<=m; j++) {

		for (k=1; k<=n; k++)  y[k] = Y[k][j];
		ldl_dcmp( A, n, diag, y, x, 0, 1, &ok ); /*  L D L'  bksbtn */

		if (disp)
		 fprintf(stderr,"    LDL' RMS matrix precision:");
		error = ok = 1;
		do {					/*improve the solution*/
			ldl_mprove ( A, n, diag, y, x, &error, &ok );
			if (disp) fprintf(stderr,"%9.2e", error );
		} while ( ok );
		if (disp) fprintf(stderr,"\n");

		for (i=1; i<=m; i++) {
			Ac[i][j] = 0.0;
			for (k=1; k<=n; k++) Ac[i][j] += X[k][i] * x[k];
		}
	}

	free_vector(diag,1,n);
	free_vector(x,1,n);
	free_vector(y,1,n);
	return;
}


/*------------------------------------------------------------------------------
AIXAI  -  calculate quadratic form with inverse matrix    inv(A) * X * inv(A)  
          A is n by n    X is n by n                                    15sep01
------------------------------------------------------------------------------*/
void aixai ( A, X, n )
float	**A, **X;
int	n;
{
	float	*diag, *b, *x, **Ai, **XAi, Aij, error;
	int	i,j,k, ok, disp=0;

	diag = vector(1,n);
	x    = vector(1,n);
	b    = vector(1,n);
	Ai   = matrix(1,n,1,n);
	XAi  = matrix(1,n,1,n);

	for (i=1; i<=n; i++) {
		diag[i] = x[i] = b[i] = 0.0;
		for (j=1; j<=n; j++)	XAi[i][j] = Ai[i][j] = 0.0;
	}

	ldl_dcmp ( A, n, diag, b, x, 1, 0, &ok );	/*  L D L'  decomp */

	for (j=1; j<=n; j++) {				/* compute inv(A)  */

		for (k=1; k<=n; k++)  b[k] = 0.0;
		b[j] = 1.0;
		ldl_dcmp( A, n, diag, b, x, 0, 1, &ok ); /*  L D L'  bksbtn */

		if (disp)
		 fprintf(stderr,"    LDL' RMS matrix precision:");
		error = ok = 1;
		do {					/*improve the solution*/
			ldl_mprove ( A, n, diag, b, x, &error, &ok );
			if (disp) fprintf(stderr,"%9.2e", error );
		} while ( ok );
		if (disp) fprintf(stderr,"\n");

		for (k=1; k<=n; k++)  Ai[j][k] = x[k];	/* save inv(A) */
	}

	for (i=1; i<=n; i++)				/* make symmetric */
		for (j=i; j<=n; j++)
			Ai[i][j] = Ai[j][i] = 0.5 * ( Ai[i][j] + Ai[j][i] );

	for (i=1; i<=n; i++) { 			/* compute X * inv(A)	*/
		for (j=1; j<=n; j++) {
			Aij = 0.0;
			for (k=1; k<=n; k++)	Aij += X[i][k]*Ai[k][j];
			XAi[i][j] = Aij;
		}
	}

	for (i=1; i<=n; i++) {		/* compute inv(A) * X * inv(A)	*/
		for (j=1; j<=n; j++) {
			Aij = 0.0;
			for (k=1; k<=n; k++)	Aij += Ai[i][k] * XAi[k][j];
			X[i][j] = Aij;
		}
	}
	for (i=1; i<=n; i++)				/* make symmetric */
		for (j=i; j<=n; j++)
			X[i][j] = X[j][i] = 0.5 * ( X[i][j] + X[j][i] );

	free_vector ( diag, 1,n );
	free_vector ( x, 1,n );
	free_vector ( b, 1,n );
	free_matrix ( Ai, 1,n,1,n );
	free_matrix ( XAi, 1,n,1,n );

	return;
}

/*------------------------------------------------------------------------------
DEALLOCATE  -  release the allocated memory				23feb94
------------------------------------------------------------------------------*/
void deallocate( x,y,z,r, L,Le, J1, J2, Ax,Asy,Asz, J,Iy,Iz, E, G,
      K,Q,F,D,R,W,P,T, feF,Fo, d,BMs,JMs,JMx,JMy,JMz, M,f,V, nJ,nM,DoF, modes ) 

int	nJ, nM, DoF, *J1, *J2, *R;
float	*x,*y,*z, *r, *L, *Le, *Ax, *Asy,*Asz, *J,*Iy,*Iz, *E, *G,  **K,
**Q,*F,*D,**W,**P,**T, **feF, *Fo, *d,*BMs,*JMs,*JMx,*JMy,*JMz, **M,*f,**V;
{
	free_vector(x,1,nJ);
	free_vector(y,1,nJ);
	free_vector(z,1,nJ);
	free_vector(r,1,nJ);
	free_vector(L,1,nM);
	free_vector(Le,1,nM);

	free_ivector(J1,1,nM);
	free_ivector(J2,1,nM);

	free_vector(Ax,1,nM);
	free_vector(Asy,1,nM);
	free_vector(Asz,1,nM);
	free_vector(J,1,nM);
	free_vector(Iy,1,nM);
	free_vector(Iz,1,nM);
	free_vector(E,1,nM);
	free_vector(G,1,nM);

	free_matrix(K,1,DoF,1,DoF);
	free_matrix(Q,1,nM,1,12);

	free_vector(F,1,DoF);
	free_vector(D,1,DoF);
	free_ivector(R,1,DoF);
	free_matrix(W,1,nM,1,4);
	free_matrix(P,1,nM,1,5);
	free_matrix(T,1,nM,1,8);
	free_matrix(feF,1,nM,1,12);
	free_vector(Fo,1,DoF);

	free_vector(d,1,nM);
	free_vector(BMs,1,nM);
	free_vector(JMs,1,nJ);
	free_vector(JMx,1,nJ);
	free_vector(JMy,1,nJ);
	free_vector(JMz,1,nJ);

	if ( modes > 0 ) {
		free_matrix(M,1,DoF,1,DoF);
		free_vector(f,1,modes);
		free_matrix(V,1,DoF,1,DoF);
	}

	return;
}

/* itoa moved to frm_io.c */

/* removed strcat -- it's in <string.h> in the standard C library */

/* removed strcpy -- it's in <string.h> in the standard C library */

/* dots moved to frm_io.c */

