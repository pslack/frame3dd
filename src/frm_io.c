/*******************************************************************************
 File frm_io.c
 File input/output routines for the '.frm' format.

 Part of FRAME -- see frame.c for details.
 ---------------------------------------------------------------------------
 Copyright (C) 1992-2007  Henri P. Gavin
 
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

 INPUT FORMAT is documented in frame.c

******************************************************************************/

#include "frm_io.h"

#include <math.h>
#include <time.h>

/*------------------------------------------------------------------------------
READ_INPUT  -  read material and geometry data, calc lengths		15dec97
------------------------------------------------------------------------------*/
void read_input(
		FILE *fp
		, int nJ, int nM, float *x, float *y, float *z
		, float *r, float *L, float *Le
		, int *J1, int *J2, int *anlyz, int *geom, float **Q
		, float *Ax, float *Asy, float *Asz
		, float *J, float *Iy, float *Iz, float *E, float *G, float *p
		, int *shear, char meshfile[], char plotfile[], float *exagg
){

	int	j1, j2, i, j, X=0, Y=0, Z=0;
	void	exit();

	for (i=1;i<=nJ;i++) {		/* read joint coordinates	*/
		fscanf(fp, "%d", &j );
		if ( j <= 0 || j > nJ ) {
		    fprintf(stderr,"  error in joint coordinate data: joint number out of range  ");
		    fprintf(stderr,"  Joint: %d  \n", j);
		    exit(1);
		}
		fscanf(fp, "%lf %lf %lf %lf", &x[j], &y[j], &z[j], &r[j]);
		r[j] = fabs(r[j]); 
	}
	for (i=1;i<=nM;i++) {		/* read member properties	*/
		fscanf(fp, "%d", &j );
		if ( j <= 0 || j > nM ) {
		    fprintf(stderr,"  error in member property data: Member number out of range  ");
		    fprintf(stderr,"  Member: %d  \n", j);
		    exit(1);
		}
		fscanf(fp, "%d %d", &J1[j], &J2[j] );
		if ( J1[j] <= 0 || J1[j] > nJ || J2[j] <= 0 || J2[j] > nJ ) {
		    fprintf(stderr,"  error in member property data: joint number out of range  ");
		    fprintf(stderr,"  Member: %d \n", j);
		    exit(1);
		}
		fscanf(fp, "%lf %lf %lf", &Ax[j], &Asy[j], &Asz[j] );
		fscanf(fp, "%lf %lf %lf", &J[j],  &Iy[j],  &Iz[j] );
		fscanf(fp, "%lf %lf %lf", &E[j], &G[j], &p[j]);

		p[j] = p[j]*PI/180.0;	/* convert from degrees to radians */

		if ( Ax[j] < 0 || Asy[j] < 0 || Asz[j] < 0 ||
		      J[j] < 0 ||  Iy[j] < 0 ||  Iz[j] < 0	) {
		    fprintf(stderr,"  error in member property data: member section property < 0  ");
		    fprintf(stderr,"  Member: %d  \n", j);
		    exit(1);
		}
		if ( Ax[j] == 0 ) {
		    fprintf(stderr,"  error in member property data: cross section area is zero   ");
		    fprintf(stderr,"  Member: %d  \n", j);
		    exit(1);
		}
		if ( (Asy[j] == 0 || Asz[j] == 0) && G[j] == 0 ) {
		    fprintf(stderr,"  error in member property data: a shear area and shear modulus are zero   ");
		    fprintf(stderr,"  Member: %d  \n", j);
		    exit(1);
		}
		if ( J[j] == 0 ) {
		    fprintf(stderr,"  error in member property data: torsional moment of inertia is zero   ");
		    fprintf(stderr,"  Member: %d  \n", j);
		    exit(1);
		}
		if ( Iy[j] == 0 || Iz[j] == 0 ) {
		    fprintf(stderr,"  error: cross section bending moment of inertia is zero   ");
		    fprintf(stderr,"  Member: %d  \n", j);
		    exit(1);
		}
		if ( E[j] <= 0 || G[j] <= 0 ) {
		    fprintf(stderr,"  error : material elastic modulus E or G f is not positive   ");
		    fprintf(stderr,"  Member: %d  \n", j);
		    exit(1);
		}
	}
	for (i=1;i<=nM;i++) {		/* calculate member lengths	*/
		j1 = J1[i];
		j2 = J2[i];
		L[i] =	(x[j2]-x[j1]) * (x[j2]-x[j1]) +
			(y[j2]-y[j1]) * (y[j2]-y[j1]) +
			(z[j2]-z[j1]) * (z[j2]-z[j1]);
		L[i] = sqrt( L[i] );
		Le[i] = L[i] - r[j1] - r[j2];
		if ( j1 == j2 || L[i] == 0.0 ) {
		   fprintf(stderr,
			" Members must start and stop at different joints\n");
		   fprintf(stderr,
			" member %d  J1= %d J2= %d L= %e\n", i, j1,j2, L[i] );
		   fprintf(stderr,
			" Perhaps member %d has not been specified. \n", i );
		   exit(1);
		}
		if ( Le[i] <= 0.0 ) {
		   fprintf(stderr, " Joint radii are too large.\n");
		   fprintf(stderr,
			" member %d  J1= %d J2= %d L= %e \n", i, j1,j2, L[i] );
		   fprintf(stderr,
			" r1= %e r2= %e Le= %e \n", r[j1], r[j2], Le[i] );
		   exit(1);
		}
	}
	fscanf( fp, "%d %d %s %s %lf %d",
			shear, geom, meshfile, plotfile, exagg, anlyz );

	if (*shear != 0 && *shear != 1) {
	    fprintf(stderr," Rember to specify shear deformations");
	    fprintf(stderr," with a 0 or a 1 after the member info.\n");
	    exit(1);
	}

	if (*geom != 0 && *geom != 1) {
	    fprintf(stderr," Rember to specify geometric stiffness");
	    fprintf(stderr," with a 0 or a 1 after the member info.\n");
	    exit(1);
	}

	if ( *exagg < 0 ) {
	    fprintf(stderr," Remember to specify an exageration");
	    fprintf(stderr," factor greater than zero\n");
	    exit(1);
	}

	for (i=1;i<=nM;i++)	for(j=1;j<=12;j++)	Q[i][j] = 0.0;

	return;
}


/*-----------------------------------------------------------------------------
GETLINE  -  get line into a character string. from K&R                  3feb94
-----------------------------------------------------------------------------*/
void getline (
		FILE	*fp
		, char    *s
		, int     lim
){
    int     c=0, i=0;

    while (--lim > 0 && (c=getc(fp)) != EOF && c != '\n' )
            s[i++] = c;
/*      if (c == '\n')  s[i++] = c;	*/
    s[i] = '\0';
    return;
}


/*-----------------------------------------------------------------------------
PARSE_INPUT                                                             7may03
 remove comments from the input file, and write a 'clean' input file 
-----------------------------------------------------------------------------*/
void parse_input(FILE *fp){
	FILE	*fpc;		/* cleaned inout/output file pointer	*/
	void	getline_no_comment();
	char	line[256];
	void	exit();

	if ((fpc = fopen ("frame.cln", "w")) == NULL) {	
		fprintf (stderr," error: cannot open file 'frame.cln'\n");
		exit(1);
	}

	do {
		getline_no_comment(fp, line, 256);
		fprintf(fpc, "%s \n", line );
	} while ( line[0] != '_' && line[0] != EOF );

	fclose(fpc);

}


/*-----------------------------------------------------------------------------
GETLINE_NO_COMMENT                                                      7may03
 get a line into a character string. from K&R
 get the line only up to one of the following characters:  \n  %  #  ;  ? 

-----------------------------------------------------------------------------*/
void getline_no_comment(
		FILE *fp            /**< pointer to the file from which to read */
		,char *s             /**< pointer to the string to which to write */
		,int lim            /**< the longest anticipated line length  */
){
    int     c=0, i=0;

    while (--lim > 0 && (c=getc(fp)) != EOF && c != '\n'
            && c != '%' && c != '#' && c != ';' && c != '?' )
            s[i++] = c;
/*      if (c == '\n')  s[i++] = c;     */
    s[i] = '\0';
	if (c != '\n') 
		while (--lim > 0 && (c=getc(fp)) != EOF && c != '\n') 
		/* read the rest of the line, otherwise do nothing */ ;

	if ( c == EOF ) s[0] = EOF;

    return;
}


/*------------------------------------------------------------------------------
READ_LOADS  -  read load information data, form un-restrained load vector	6dec06
------------------------------------------------------------------------------*/
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
){
	float	Nx1, Vy1, Vz1, Mx1, My1, Mz1,	/* fixed end forces */
		Nx2, Vy2, Vz2, Mx2, My2, Mz2,
		Ksy, Ksz, 		/* shear deformatn coefficients	*/
		a, b,				/* point load locations */
		hy, hz,			/* section dimensions in local coords */
		t1, t2, t3, t4, t5, t6, t7, t8, t9;	/* 3D coord Xfrm coef */
	int	i,j,l,n, j1, j2;
	void dots(),
		exit();

	for (j=1; j<=DoF; j++)
		F_mech[j] = F_temp[j] = 0.0;
	for (i=1; i<=12; i++) 
		for (n=1; n<=nM; n++) 
			feF_mech[n][i] = feF_temp[n][i] = 0.0;

	fscanf(fp,"%d", nF );		/* joint point loads		*/
	printf(" number of loaded joints ");
	dots(28);
	printf(" nF = %d\n",*nF);
	for (i=1; i <= *nF; i++) {	/* ! global structural coordinates ! */
		fscanf(fp,"%d", &j);
		if ( j < 1 || j > nJ ) {
		    fprintf(stderr,"  error in joint load data: joint number out of range  ");
		    fprintf(stderr,"  Joint: %d  \n", j);
		    fprintf(stderr,"  Perhaps you did not specify %d joint loads \n", *nF );
		    exit(1);
		}
		for (l=5; l>=0; l--)	fscanf(fp,"%lf", &F_mech[6*j-l] );
		if ( F_mech[6*j-5]==0 && F_mech[6*j-4]==0 && F_mech[6*j-3]==0 && F_mech[6*j-2]==0 && F_mech[6*j-1] && F_mech[6*j]==0 )
		    fprintf(stderr,"  warning: All joint loads applied at joint %d  are zero\n", j );
	}

	fscanf(fp,"%d", nW );		/* uniform distributed loads	*/
	printf(" number of uniform distributed loads ");
	dots(16);
	printf(" nW = %d\n",*nW);
	if ( *nW < 0 || *nW > nM ) {
		fprintf(stderr,"  error: valid ranges for nW is 0 ... %d \n", nM );
		exit(1);
	}
	for (i=1; i <= *nW; i++) {	/* ! local member coordinates !	*/
		fscanf(fp,"%d", &n );		
		if ( n < 1 || n > nM ) {
		    fprintf(stderr,"  error in uniform distributed loads: member number %d is out of range\n",n);
		    exit(1);
		}
		W[i][1] = (float) n;
		for (l=2; l<=4; l++)	fscanf(fp,"%lf", &W[i][l] );

		if ( W[i][2]==0 && W[i][3]==0 && W[i][4]==0 )
		    fprintf(stderr,"  warning: All distributed loads applied to member %d  are zero\n", (int)W[i][1] );

		Nx1 = Nx2 = W[i][2]*Le[n] / 2.0;
		Vy1 = Vy2 = W[i][3]*Le[n] / 2.0;
		Vz1 = Vz2 = W[i][4]*Le[n] / 2.0;
		Mx1 = Mx2 = 0.0;
		My1 = -W[i][4]*Le[n]*Le[n] / 12.0;	My2 = -My1;
		Mz1 =  W[i][3]*Le[n]*Le[n] / 12.0;	Mz2 = -Mz1;

		/* debugging
		printf("n=%d Vy=%9.2e Vz=%9.2e My=%9.2e Mz=%9.2e\n",
						n, Vy1,Vz1, My1,Mz1 );	*/

		j1 = J1[n];	j2 = J2[n];

		coord_trans ( x, y, z, L[n], j1, j2,
			&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p[n] );

		/* debugging
		printf("t1=%5.2f t2=%5.2f t3=%5.2f \n", t1, t2, t3 );
                printf("t4=%5.2f t5=%5.2f t6=%5.2f \n", t4, t5, t6 );
                printf("t7=%5.2f t8=%5.2f t9=%5.2f \n", t7, t8, t9 ); */

		/* {F} = [T]'{Q} */
		feF_mech[n][1]  += ( Nx1*t1 + Vy1*t4 + Vz1*t7 );
		feF_mech[n][2]  += ( Nx1*t2 + Vy1*t5 + Vz1*t8 );
		feF_mech[n][3]  += ( Nx1*t3 + Vy1*t6 + Vz1*t9 );
		feF_mech[n][4]  += ( Mx1*t1 + My1*t4 + Mz1*t7 );
		feF_mech[n][5]  += ( Mx1*t2 + My1*t5 + Mz1*t8 );
		feF_mech[n][6]  += ( Mx1*t3 + My1*t6 + Mz1*t9 );

		feF_mech[n][7]  += ( Nx2*t1 + Vy2*t4 + Vz2*t7 );
		feF_mech[n][8]  += ( Nx2*t2 + Vy2*t5 + Vz2*t8 );
		feF_mech[n][9]  += ( Nx2*t3 + Vy2*t6 + Vz2*t9 );
		feF_mech[n][10] += ( Mx2*t1 + My2*t4 + Mz2*t7 );
		feF_mech[n][11] += ( Mx2*t2 + My2*t5 + Mz2*t8 );
		feF_mech[n][12] += ( Mx2*t3 + My2*t6 + Mz2*t9 );

		/* debugging
		printf("n=%d ", n);
		for (l=1;l<=12;l++) {
			if (feF_mech[n][l] != 0)
				printf(" feF %d = %9.2e ", l, feF_mech[n][l] );
		}
		printf("\n"); */
	}

	fscanf(fp,"%d", nP );		/* concentrated point loads	*/
	printf(" number of concentrated member point loads ");
	dots(10);
	printf(" nP = %d\n",*nP);
	if ( *nP < 0 || *nP > nM ) {
		fprintf(stderr,"  error: valid ranges for nP is 0 ... %d \n", nM );
		exit(1);
	}
	for (i=1; i <= *nP; i++) {	/* ! local member coordinates !	*/
		fscanf(fp,"%d", &n );		
		if ( n < 1 || n > nM ) {
		    fprintf(stderr,"  error in internal point loads: member number %d is out of range\n",n);
		    exit(1);
		}
		P[i][1] = (float) n;
		for (l=2; l<=5; l++)	fscanf(fp,"%lf", &P[i][l] );
		a = P[i][5];	b = L[n] - a;

		if ( a < 0 || L[n] < a || b < 0 || L[n] < b ) {
		    fprintf(stderr,"  error in point load data: Point load coord. out of range\n");
		    fprintf(stderr,"  Member: %d  L: %lf  load coord.: %lf\n",
							n, L[n], P[i][5] );
		    exit(1);
		}

		if ( shear ) {
			Ksy = G[n]*Asy[n]*Le[n]*Le[n] / (12.*E[n]*Iz[n]);
			Ksz = G[n]*Asz[n]*Le[n]*Le[n] / (12.*E[n]*Iy[n]);
		} else	Ksy = Ksz = 0.0;


		Nx1 = P[i][2]*a/L[n];
		Nx2 = P[i][2]*b/L[n];
		Vy1 = (1./(1.+Ksz))*P[i][3]*b*b*(3.*a + b) / ( L[n]*L[n]*L[n] )+
			(Ksz/(1.+Ksz)) * P[i][3]*b/L[n];
		Vy2 = (1./(1.+Ksz))*P[i][3]*a*a*(3.*b + a) / ( L[n]*L[n]*L[n] )+
			(Ksz/(1.+Ksz)) * P[i][3]*a/L[n];
		Vz1 = (1./(1.+Ksy))*P[i][4]*b*b*(3.*a + b) / ( L[n]*L[n]*L[n] )+
			(Ksy/(1.+Ksy)) * P[i][4]*b/L[n];
		Vz2 = (1./(1.+Ksy))*P[i][4]*a*a*(3.*b + a) / ( L[n]*L[n]*L[n] )+
			(Ksy/(1.+Ksy)) * P[i][4]*a/L[n];
		Mx1 = Mx2 = 0.0;
		My1 = -(1./(1.+Ksy)) * P[i][4]*a*b*b / ( L[n]*L[n] ) - 
			(Ksy/(1.+Ksy))* P[i][4]*a*b / (2.*L[n]);
		My2 =  (1./(1.+Ksy)) * P[i][4]*a*a*b / ( L[n]*L[n] ) + 
			(Ksy/(1.+Ksy))* P[i][4]*a*b / (2.*L[n]);
		Mz1 =  (1./(1.+Ksz)) * P[i][3]*a*b*b / ( L[n]*L[n] ) + 
			(Ksz/(1.+Ksz))* P[i][3]*a*b / (2.*L[n]);
		Mz2 = -(1./(1.+Ksz)) * P[i][3]*a*a*b / ( L[n]*L[n] ) - 
			(Ksz/(1.+Ksz))* P[i][3]*a*b / (2.*L[n]);

		j1 = J1[n];	j2 = J2[n];

		coord_trans ( x, y, z, L[n], j1, j2,
			&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p[n] );

		/* {F} = [T]'{Q} */
		feF_mech[n][1]  += ( Nx1*t1 + Vy1*t4 + Vz1*t7 );
		feF_mech[n][2]  += ( Nx1*t2 + Vy1*t5 + Vz1*t8 );
		feF_mech[n][3]  += ( Nx1*t3 + Vy1*t6 + Vz1*t9 );
		feF_mech[n][4]  += ( Mx1*t1 + My1*t4 + Mz1*t7 );
		feF_mech[n][5]  += ( Mx1*t2 + My1*t5 + Mz1*t8 );
		feF_mech[n][6]  += ( Mx1*t3 + My1*t6 + Mz1*t9 );

		feF_mech[n][7]  += ( Nx2*t1 + Vy2*t4 + Vz2*t7 );
		feF_mech[n][8]  += ( Nx2*t2 + Vy2*t5 + Vz2*t8 );
		feF_mech[n][9]  += ( Nx2*t3 + Vy2*t6 + Vz2*t9 );
		feF_mech[n][10] += ( Mx2*t1 + My2*t4 + Mz2*t7 );
		feF_mech[n][11] += ( Mx2*t2 + My2*t5 + Mz2*t8 );
		feF_mech[n][12] += ( Mx2*t3 + My2*t6 + Mz2*t9 );
	}

	fscanf(fp,"%d", nT );		/* thermal loads		*/
	printf(" number of members with temperature changes ");
	dots(9);
	printf(" nT = %d\n",*nT);
	if ( *nT < 0 || *nT > nM ) {
		fprintf(stderr,"  error: valid ranges for nT is 0 ... %d \n", nM );
		exit(1);
	}
	for (i=1; i <= *nT; i++) {	/* ! member coordinates !	*/
		fscanf(fp,"%d", &n );		
		if ( n < 1 || n > nM ) {
		    fprintf(stderr,"  error in temperature loads: member number %d is out of range\n",n);
		    exit(1);
		}
		T[i][1] = (float) n;
		for (l=2; l<=8; l++)	fscanf(fp,"%lf", &T[i][l] );
		a  = T[i][2];
		hy = T[i][3];
		hz = T[i][4];

		if ( hy < 0 || hz < 0 ) {
		    fprintf(stderr,"  error in thermal load data: section dimension < 0\n");
		    fprintf(stderr,"  Member: %d  hy: %lf  hz: %lf\n", n,hy,hz);
		    exit(1);
		}

		Nx2 = (a/4.0)*( T[i][5]+T[i][6]+T[i][7]+T[i][8])*E[n]*Ax[n];
		Nx1 = -Nx2;
		Vy1 = Vy2 = Vz1 = Vz2 = 0.0; 
		Mx1 = Mx2 = 0.0;
		My1 =  (a/hz)*(T[i][8]-T[i][7])*E[n]*Iy[n];
		My2 = -My1;
		Mz1 =  (a/hy)*(T[i][5]-T[i][6])*E[n]*Iz[n];
		Mz2 = -Mz1;

		j1 = J1[n];	j2 = J2[n];

		coord_trans ( x, y, z, L[n], j1, j2,
			&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p[n] );

		/* {F} = [T]'{Q} */
		feF_temp[n][1]  += ( Nx1*t1 + Vy1*t4 + Vz1*t7 );
		feF_temp[n][2]  += ( Nx1*t2 + Vy1*t5 + Vz1*t8 );
		feF_temp[n][3]  += ( Nx1*t3 + Vy1*t6 + Vz1*t9 );
		feF_temp[n][4]  += ( Mx1*t1 + My1*t4 + Mz1*t7 );
		feF_temp[n][5]  += ( Mx1*t2 + My1*t5 + Mz1*t8 );
		feF_temp[n][6]  += ( Mx1*t3 + My1*t6 + Mz1*t9 );

		feF_temp[n][7]  += ( Nx2*t1 + Vy2*t4 + Vz2*t7 );
		feF_temp[n][8]  += ( Nx2*t2 + Vy2*t5 + Vz2*t8 );
		feF_temp[n][9]  += ( Nx2*t3 + Vy2*t6 + Vz2*t9 );
		feF_temp[n][10] += ( Mx2*t1 + My2*t4 + Mz2*t7 );
		feF_temp[n][11] += ( Mx2*t2 + My2*t5 + Mz2*t8 );
		feF_temp[n][12] += ( Mx2*t3 + My2*t6 + Mz2*t9 );
	}

	for (n=1; n<=nM; n++) {
		j1 = J1[n];	j2 = J2[n];
		for (i=1; i<= 6; i++)	F_mech[6*j1- 6+i] += feF_mech[n][i];
		for (i=7; i<=12; i++)	F_mech[6*j2-12+i] += feF_mech[n][i];
		for (i=1; i<= 6; i++)	F_temp[6*j1- 6+i] += feF_temp[n][i];
		for (i=7; i<=12; i++)	F_temp[6*j2-12+i] += feF_temp[n][i];
	}

	return;
}


/*------------------------------------------------------------------------------
READ_REACTIONS  -  Read fixed joint displacement boundary conditions	27dec01
------------------------------------------------------------------------------*/
void read_reactions (
		FILE *fp
		, int DoF, int *nD, int *nR
		, int nJ, float *Dp, int *R, int *sumR 
){
	int	i,j,l;
	void	dots(), exit();

	for (i=1; i<=DoF; i++)  {
		Dp[i] = 0.0;
		R[i] = 0;
	}

	fscanf(fp,"%d", nR );	/* read restrained degrees of freedom */
	printf(" number of joints with reactions ");
	dots(20);
	printf(" nR = %d\n",*nR);
	if ( *nR < 0 || *nR > DoF/6 ) {
		fprintf(stderr,"  error: valid ranges for nR is 0 ... %d \n", DoF/6 );
		exit(1);
	}

	for (i=1; i <= *nR; i++) {
	    fscanf(fp,"%d", &j);
	    for (l=5; l >=0; l--) {

		fscanf(fp,"%d", &R[6*j-l] );

		if ( j > nJ ) {
		    fprintf(stderr,"  error in reaction data: joint number %d is greater than the number of joints, %d \n", j, nJ );
		    exit(1);
		}
		if ( R[6*j-l] != 0 && R[6*j-l] != 1 ) {
		    fprintf(stderr,"  error in reaction data: Reaction data must be 0 or 1\n");
		    fprintf(stderr,"  Data for joint %d, DoF %d is %d\n",
							j, 6-l, R[6*j-l] );
		    exit(1);
		}
	    }
	    *sumR = 0;
	    for (l=5; l >=0; l--) 	*sumR += R[6*j-l];
	    if ( *sumR == 0 ) {
		fprintf(stderr,"  error: joint %3d has no reactions\n", j);
		fprintf(stderr,"  Remove joint %3d from the list of reactions\n", j);
		fprintf(stderr,"  and set nR to %3d \n", *nR-1);
		exit(1);
	    }
	}
	*sumR=0;	for (i=1;i<=DoF;i++)	*sumR += R[i];
	if ( *sumR < 4 ) {
	    fprintf(stderr,"  warning:  Un-restrained structure\n");
	    fprintf(stderr,"  %d imposed reactions.", *sumR ); 
	    fprintf(stderr,"  At least 4 reactions are required.\n");
	    /*	exit(1); */
	}
	if ( *sumR >= DoF ) {
	    fprintf(stderr,"  error in reaction data:  Fully restrained structure\n");
	    fprintf(stderr,"  %d imposed reactions >= %d degrees of freedom\n",
								*sumR, DoF );
	    exit(1);
	}

	fscanf(fp,"%d", nD );		/* read prescribed displacements */
	printf(" number of joints with prescribed displacements ");
	dots(5);
	printf(" nD = %d\n",*nD);
	for (i=1; i <= *nD; i++) {
		fscanf(fp,"%d", &j);
		for (l=5; l >=0; l--) {
			fscanf(fp,"%lf", &Dp[6*j-l] );
			if ( R[6*j-l] == 0 && Dp[6*j-l] != 0.0 ) {
			    printf(" Initial displacements can be prescribed");
			    printf(" only at restrained coordinates\n");
			    printf(" joint: %d  dof: %d  R: %d\n",
							j, 6-l, R[6*j-l] );
			    exit(1);
			}
		}
	}
	return;
}

/*------------------------------------------------------------------------------
READ_MASSES  -  read member densities and extra inertial mass data	16aug01
------------------------------------------------------------------------------*/
void read_masses(
		FILE *fp
		, int nJ, int nM, int *nI
		, float *d, float *BMs, float *JMs, float *JMx, float *JMy, float *JMz
		, float *L, float *Ax
		, float *total_mass, float *struct_mass
		, int *modes, int *Mmethod, int *lump
		, char modefile[]
		, float *tol, float *shift, int anim[], int *pan
){
	FILE	*mf;				/* mass data file	*/
	float	ms = 0.0;
	int	chk, j, jnt, m, mem, nA;
	void	dots(), exit();

	*total_mass = *struct_mass = 0.0;	


	chk = fscanf ( fp, "%d", modes );

	printf(" number of dynamic modes ");
	dots(28);
	printf(" modes = %d\n", *modes);

	if ( *modes < 1 || chk != 1 ) {
		*modes = 0;
		return;
	}

	fscanf( fp, "%d", Mmethod );

	printf(" modal analysis method ");
	dots(30);
	printf(" %d ",*Mmethod);
	if ( *Mmethod == 1 ) printf(" (Subspace-Jacobi)\n");
	if ( *Mmethod == 2 ) printf(" (Stodola)\n");


	/*
	mf = fopen("MassData.txt","w");	// open mass data file 
	if ((mf = fopen ("MassData.txt", "w")) == NULL) {	
	  fprintf (stderr," error: cannot open file 'MassData.txt'\n");
	  exit(1);
	}
	fprintf(mf,"%% structural mass data \n");
	fprintf(mf,"%% element\tAx\t\tlength\t\tdensity\t\tmass \n");
	*/

	fscanf( fp, "%d", lump );
	fscanf( fp, "%s", modefile );
	fscanf( fp, "%lf", tol );
	fscanf( fp, "%lf", shift );
	for (m=1; m <= nM; m++) {	/* read inertia data	*/
		fscanf(fp, "%d", &mem );
		fscanf(fp, "%lf %lf", &d[mem], &BMs[mem] );
		*total_mass  += d[mem]*Ax[mem]*L[mem] + BMs[mem];
		*struct_mass += d[mem]*Ax[mem]*L[mem];
		/*
		fprintf(mf," %4d\t\t%12.5e\t%12.5e\t%12.5e\t%12.5e \n",
		 mem, Ax[mem], L[mem], d[mem], d[mem]*Ax[mem]*L[mem] );
		*/
	}

	/*
	fclose(mf);
	*/

	/* number of joints with extra inertias */
	fscanf(fp,"%d", nI );
	printf(" number of joints with extra lumped inertia ");
        dots(9);
        printf(" nI = %d\n",*nI);
	for (j=1; j <= *nI; j++) {
		fscanf(fp, "%d", &jnt );
		if ( jnt < 1 || jnt > nJ ) {
	    		fprintf(stderr,"  error in joint load data: joint number out of range  ");
	    		fprintf(stderr,"  Joint: %d  \n", j);
	    		fprintf(stderr,"  Perhaps you did not specify %d extra masses \n", *nI );
	    		exit(1);
		}
		fscanf(fp, "%lf %lf %lf %lf",
			&JMs[jnt], &JMx[jnt], &JMy[jnt], &JMz[jnt] );
		*total_mass += JMs[jnt];

		if ( JMs[jnt]==0 && JMx[jnt]==0 && JMy[jnt]==0 && JMz[jnt]==0 )
	    	fprintf(stderr,"  warning: All extra joint inertia at joint %d  are zero\n", jnt );
	}

	for (m=1;m<=nM;m++) {			/* check inertia data	*/
	    if ( d[m] < 0.0 || BMs[m] < 0.0 || d[m]+BMs[m] <= 0.0 ) {
		fprintf(stderr,"  error: Non-positive mass or density\n");
		fprintf(stderr,"  d[%d]= %lf  BMs[%d]= %lf\n",m,d[m],m,BMs[m]);
		exit(1);
	    }
	}
/*	for (m=1;m<=nM;m++)	ms += BMs[m];/* consistent mass doesn't agree */
/*	if ( ms > 0.0 )		*lump = 1;  /* with concentrated masses, BMs  */

	printf(" structural mass ");
	dots(36);
	printf("  %12.4e\n",*struct_mass);
	printf(" total mass ");
	dots(41);
	printf("  %12.4e\n",*total_mass);
	fscanf ( fp, "%d", &nA );
	printf(" number of modes to be animated ");
	dots(21);
	printf(" nA = %d\n",nA);
	if (nA > 20) 
	  printf(" nA = %d, only 20 or fewer modes may be animated\n", nA );
	for ( m = 0; m < 20; m++ )	anim[m] = 0;
	for ( m = 0; m < nA; m++ )	fscanf ( fp, "%d", &anim[m] );

	fscanf ( fp, "%d", pan );

	return;
}

/*------------------------------------------------------------------------------
READ_CONDENSE   -  read matrix condensation information 	        30aug01
------------------------------------------------------------------------------*/
void read_condense (
		FILE *fp
		, int nJ, int modes
		, int *nC, int *Cdof, int *Cmethod, int *q, int *m
){
	int	i,j,k,  chk, **qm, **imatrix();
	void	dots(),
		free_imatrix(), exit();

	*Cmethod = *nC = *Cdof = 0;

	if ( (chk = fscanf ( fp, "%d", Cmethod )) != 1 )   {
		*Cmethod = *nC = *Cdof = 0;
		return;
	}

	if ( *Cmethod <= 0 )  {
		*Cmethod = *nC = *Cdof = 0;
		return;
	}

	if ( *Cmethod > 3 ) *Cmethod = 1;	/* default */
	printf(" condensation method ");
	dots(32);
	printf(" %d ", *Cmethod );
	if ( *Cmethod == 1 )	printf(" (static only) \n");
	if ( *Cmethod == 2 )	printf(" (Guyan) \n");
	if ( *Cmethod == 3 )	printf(" (dynamic) \n");

	if ( (chk = fscanf ( fp, "%d", nC )) != 1 )  {
		*Cmethod = *nC = *Cdof = 0;
		return;
	}

	printf(" number of joints with condensed DoF's ");
	dots(14);
	printf(" nC = %d\n", *nC );

	if ( (*nC) > nJ ) {
	  fprintf(stderr," error in matrix condensation data: \n");
	  fprintf(stderr,"  error: nC > nJ ... nC=%d; nJ=%d;\n",*nC,nJ);
	  fprintf(stderr,"  The number of joints with condensed DoF's ");
	  fprintf(stderr,"may not exceed the total number of joints.\n");
	  exit(1);
	}


	qm = imatrix( 1, *nC, 1,7 );

	for ( i=1; i <= *nC; i++) {
	 fscanf( fp, "%d %d %d %d %d %d %d",
	 &qm[i][1],
	 &qm[i][2], &qm[i][3], &qm[i][4], &qm[i][5], &qm[i][6], &qm[i][7]);
	 if ( qm[i][1] < 1 || qm[i][1] > nJ ) {		/* error check */
	  fprintf(stderr," error in matrix condensation data: ");
	  fprintf(stderr," condensed joint number out of range\n");
	  fprintf(stderr,"  cj[%d] = %d  ... nJ = %d  \n", i, qm[i][1], nJ );
	  exit(1);
	 }
	}

	for (i=1; i <= *nC; i++)  for (j=2; j<=7; j++)  if (qm[i][j]) (*Cdof)++;

	k=1;
	for (i=1; i <= *nC; i++) {
		for (j=2; j<=7; j++) {
			if (qm[i][j]) {
				q[k] = 6*(qm[i][1]-1) + j-1;
				++k;
			}
		}
	}

	for (i=1; i<= *Cdof; i++) {
	 fscanf( fp, "%d", &m[i] );
	 if ( (m[i] < 0 || m[i] > modes) && *Cmethod == 3 ) {
	  fprintf(stderr," error in matrix condensation data: \n");
	  fprintf(stderr,"  error: m[%d] = %d \n",i,m[i]);
	  fprintf(stderr,"  The condensed mode number must be between ");
	  fprintf(stderr,"  1 and %d (modes).\n", modes);
	  exit(1);
         }
	}

	free_imatrix(qm,1, *nC, 1,7);
	return;
}


/*------------------------------------------------------------------------------
CONTROL_DATA  -  save input data					7nov02
------------------------------------------------------------------------------*/
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
){
	int	i,j,n;
        time_t  now;            /* modern time variable type    (DJGPP) */

	(void) time(&now);

	fprintf(fp,"\n");
	for (i=1; i<=80; i++)	fprintf(fp,"_");
	fprintf(fp,"\n-- FRAME version:   20 Dec 2007,");
	fprintf(fp," GPL Copyright (C) 1992-2007, Henri P. Gavin --\n");
	fprintf(fp,"                     http://www.duke.edu/~hpgavin/frame/ \n");
	fprintf(fp," FRAME is distributed in the hope that it will be useful");
	fprintf(fp," but with no warranty;\n");
	fprintf(fp," for details see the GNU Public Licence:");
	fprintf(fp," http://www.fsf.org/copyleft/gpl.html\n");
	for (i=1; i<=80; i++)	fprintf(fp,"_"); fprintf(fp,"\n\n");
	fprintf(fp,"%s\n",title);
	fprintf(fp, "%s", ctime(&now) );
	for (i=1; i<=80; i++)	fprintf(fp,"_"); fprintf(fp,"\n");
	

	fprintf(fp,"JOINTS: %d    MEMBERS: %d   FIXED JOINTS: %d", nJ,nM,nR);
	fprintf(fp,"   PRESCRIBED DISPLACEMENTS: %d\n", nD );
	fprintf(fp,"JOINT LOADS: %d   UNIFORM MEMBER LOADS: %d   ", nF,nW );
	fprintf(fp,"CONCENTRATED MEMBER LOADS: %d   \n\n", nP );
	fprintf(fp,"For 2D problems, the Y-axis is vertical. \n");
#if Zvert
	fprintf(fp,"For 3D problems, the Z-axis is vertical. \n");
#else
	fprintf(fp,"For 3D problems, the Y-axis is vertical. \n");
#endif

	for (i=1; i<=80; i++)	fprintf(fp,"_");	fprintf(fp,"\n");

	fprintf(fp,"J O I N T   D A T A     ");
	fprintf(fp,"                                    R E S T R A I N T S\n");
	fprintf(fp,"  Joint      X              Y              Z");
	fprintf(fp,"         radius  Fx Fy Fz Mx My Mz\n");
	for (i=1; i<=nJ; i++) {
	 j = 6*(i-1);
	 fprintf(fp,"%5d %14.6f %14.6f %14.6f %8.3f  %2d %2d %2d %2d %2d %2d\n",
		i, x[i], y[i], z[i], r[i], 
			R[j+1], R[j+2], R[j+3], R[j+4], R[j+5], R[j+6] );
	}
	fprintf(fp,"M E M B E R   D A T A\t\t\t\t\t\t\t(local)\n");
	fprintf(fp,"  Member J1    J2     Ax   Asy   Asz    ");
	fprintf(fp,"Jxx     Iyy     Izz       E       G roll\n");
	for (i=1; i<= nM; i++) {
		fprintf(fp,"%5d %5d %5d %6.1f %5.1f %5.1f",
					i, J1[i],J2[i], Ax[i], Asy[i], Asz[i] );
		fprintf(fp," %6.1f %7.1f %7.1f %8.1f %7.1f %3.0f\n",
				J[i], Iy[i], Iz[i], E[i], G[i], p[i]*180.0/PI );
	}
	if ( shear )	fprintf(fp,"  Include shear deformations.\n");
	else		fprintf(fp,"  Neglect shear deformations.\n");
	if ( geom )	fprintf(fp,"  Include geometric stiffness.\n");
	else		fprintf(fp,"  Neglect geometric stiffness.\n");

	if ( nF > 0 || nW > 0 || nP > 0 || nT > 0) {
	    fprintf(fp,"J O I N T   L O A D S");
	    fprintf(fp,"  +  E Q U I V A L E N T   J O I N T   L O A D S\t(global)\n");
	    fprintf(fp,"  Joint       Fx          Fy          Fz");
	    fprintf(fp,"          Mxx         Myy         Mzz\n");
	    for (j=1; j<=nJ; j++) {
		i = 6*(j-1);
		if ( F[i+1] != 0.0 || F[i+2] != 0.0 || F[i+3] != 0.0 ||
		     F[i+4] != 0.0 || F[i+5] != 0.0 || F[i+6] != 0.0 ) {
			fprintf(fp, " %5d", j);
			for (i=5; i>=0; i--) fprintf(fp, " %11.3f", F[6*j-i] );
			fprintf(fp, "\n");
		}
	    }  
	}

	if ( nW > 0 ) {
	    fprintf(fp,"U N I F O R M   M E M B E R   L O A D S");
	    fprintf(fp,"\t\t\t\t\t(local)\n");
	    fprintf(fp,"  Member      Wx               Wy               Wz\n");
	    for (n=1; n<=nW; n++) {
		fprintf(fp, " %5d", (int) (W[n][1]) );
		for (i=2; i<=4; i++) fprintf(fp, " %16.8f", W[n][i] );
		fprintf(fp, "\n");
	    }  
	}

	if ( nP > 0 ) {
	    fprintf(fp,"C O N C E T R A T E D   P O I N T   L O A D S");
	    fprintf(fp,"\t\t\t\t(local)\n");
	    fprintf(fp,"  Member      Px          Py          Pz          x\n");
	    for (n=1; n<=nP; n++) {
		fprintf(fp, " %5d", (int) (P[n][1]) );
		for (i=2; i<=5; i++) fprintf(fp, " %11.3f", P[n][i] );
		fprintf(fp, "\n");
	    }  
	}

	if ( nT > 0 ) {
	    fprintf(fp,"M E M B E R   T E M P E R A T U R E   C H A N G E S");
	    fprintf(fp,"\t\t\t(local)\n");
	    fprintf(fp,"  Member    coef      hy        hz");
	    fprintf(fp,"        Ty+       Ty-       Tz+       Tz-\n");
	    for (n=1; n<=nT; n++) {
		fprintf(fp, " %5d", (int) (T[n][1]) );
		fprintf(fp, " %9.2e", T[n][2] );
		for (i=3; i<=8; i++) fprintf(fp, " %9.3f", T[n][i] );
		fprintf(fp, "\n");
	    }  
	}

	if ( nD > 0 ) {
	    fprintf(fp,"P R E S C R I B E D   D I S P L A C E M E N T S");
	    fprintf(fp,"\t\t\t(global)\n");
	    fprintf(fp,"  Joint       Dx          Dy          Dz");
	    fprintf(fp,"          Dxx         Dyy         Dzz\n");
	    for (j=1; j<=nJ; j++) {
		i = 6*(j-1);
		if ( Dp[i+1] != 0.0 || Dp[i+2] != 0.0 || Dp[i+3] != 0.0 ||
		     Dp[i+4] != 0.0 || Dp[i+5] != 0.0 || Dp[i+6] != 0.0 ) {
			fprintf(fp, " %5d", j);
			for (i=5; i>=0; i--) fprintf(fp, " %11.3f", Dp[6*j-i] );
			fprintf(fp, "\n");
		}
	    }
	}
	if (anlyz) {
	 fprintf(fp,"\nE L A S T I C   S T I F F N E S S   A N A L Y S I S");
	 fprintf(fp,"   via  L D L'  decomposition\n\n");
	}
	else		fprintf(fp,"D A T A   C H E C K   O N L Y\n");
	fflush(fp);
	return;
}


/*------------------------------------------------------------------------------
SAVE_RESULTS -  save joint displacements and member end forces		15oct98
------------------------------------------------------------------------------*/
void save_results (
		FILE *fp, int nJ, int nM, int DoF, int *J1, int *J2
		, float *F, float *D, int *R
		, float **Q, float err
		, int ok 
){
	float	disp;
	int	i,j,n;


	if ( ok < 0 ) {
	 fprintf(fp,"  * The Stiffness Matrix is not positive-definite *\n");
	 fprintf(fp,"    Check that all six rigid-body translations are restrained\n");
	 fprintf(fp,"    If geometric stiffness is included, reduce the loads.\n");
/*	 return; */
	}
		 
	fprintf(fp,"J O I N T   D I S P L A C E M E N T S");
	fprintf(fp,"\t\t\t\t\t(global)\n");
	fprintf(fp,"  Joint    X-dsp       Y-dsp       Z-dsp");
	fprintf(fp,"       X-rot       Y-rot       Z-rot\n");
	for (j=1; j<= nJ; j++) {
	    disp = 0.0;
	    for ( i=5; i>=0; i-- ) disp += fabs( D[6*j-i] );
	    if ( disp > 0.0 ) {
		fprintf(fp," %5d", j);
		for ( i=5; i>=0; i-- ) {
                        if ( fabs(D[6*j-i]) < 1.e-8 )
                                fprintf (fp, "    0.0     ", D[6*j-i] );
                        else    fprintf (fp, " %11.6f",  D[6*j-i] );
		}
		fprintf(fp,"\n");
	    }
	}
	fprintf(fp,"M E M B E R   E N D   F O R C E S");
	fprintf(fp,"\t\t\t\t\t(local)\n");
	fprintf(fp,"  Member Joint      Nx          Vy         Vz");
	fprintf(fp,"         Txx        Myy        Mzz\n");
	for (n=1; n<= nM; n++) {
		fprintf(fp," %5d  %5d", n, J1[n]);
		if ( fabs(Q[n][1]) < 0.0001 )
			fprintf (fp, "      0.0   ");
		else    fprintf (fp, " %10.3f", Q[n][1] );
		if ( Q[n][1] >=  0.0001 ) fprintf(fp, "c");
		if ( Q[n][1] <= -0.0001 ) fprintf(fp, "t");
		for (i=2; i<=6; i++) {
			if ( fabs(Q[n][i]) < 0.0001 )
				fprintf (fp, "      0.0  ");
			else    fprintf (fp, " %10.3f", Q[n][i] );
                }
		fprintf(fp,"\n");
		fprintf(fp," %5d  %5d", n, J2[n]);
		if ( fabs(Q[n][7]) < 0.0001 )
			fprintf (fp, "      0.0   ");
		else    fprintf (fp, " %10.3f", Q[n][7] );
		if ( Q[n][7] >=  0.0001 ) fprintf(fp, "t");
		if ( Q[n][7] <= -0.0001 ) fprintf(fp, "c");
		for (i=8; i<=12; i++) {
			if ( fabs(Q[n][i]) < 0.0001 )
				fprintf (fp, "      0.0  ");
			else    fprintf (fp, " %10.3f", Q[n][i] );
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"R E A C T I O N S\t\t\t\t\t\t\t(global)\n");
	fprintf(fp,"  Joint       Fx          Fy          Fz");
	fprintf(fp,"         Mxx         Myy         Mzz\n");
	for (j=1; j<=nJ; j++) {
		i = 6*(j-1);
		if ( R[i+1] || R[i+2] || R[i+3] ||
		     R[i+4] || R[i+5] || R[i+6] ) {
			fprintf(fp, " %5d", j);
                        for (i=5; i>=0; i--) {
                                if ( !R[6*j-i] || fabs(F[6*j-i]) < 0.0001 )
                                        fprintf (fp, "       0.0  ");
                                else    fprintf (fp, " %11.3f", -F[6*j-i] );
                        }
			fprintf(fp, "\n");
		}
	}
	fprintf(fp,"R M S   E Q U I L I B R I U M    E R R O R: %9.3e\n", err );
	fflush(fp);
	return;
}


/*------------------------------------------------------------------------------
MODAL_RESULTS -  save modal frequencies and mode shapes			16aug01
------------------------------------------------------------------------------*/
void modal_results(
		FILE *fp
		, int nJ, int nM, int nI, int DoF
		, float **M, float *f, float **V
		, float total_mass, float struct_mass
		, int iter, int sumR, int modes
		, float shift, int lump, float tol
		, int ok
){
	int	i, j, k, m, num_modes;
	float	mpfX, mpfY, mpfZ,	/* mode participation factors	*/
		*msX, *msY, *msZ, *vector();
	float	fs;
	void	free_vector();

	msX = vector(1,DoF);
	msY = vector(1,DoF);
	msZ = vector(1,DoF);

	for (i=1; i<=DoF; i++) {
		msX[i] = msY[i] = msZ[i] = 0.0;
		for (j=1; j<=DoF; j+=6) msX[i] += M[i][j];
		for (j=2; j<=DoF; j+=6) msY[i] += M[i][j];
		for (j=3; j<=DoF; j+=6) msZ[i] += M[i][j];
	}
	
	if ( (DoF - sumR) > modes )	num_modes = modes;
	else	num_modes = DoF - sumR;

	fprintf(fp,"\nM O D A L   A N A L Y S I S   R E S U L T S\n");
	fprintf(fp,"  Total Mass:  %e   ", total_mass );
	fprintf(fp,"  Structural Mass:  %e \n", struct_mass );
	fprintf(fp,"J O I N T   M A S S E S");
	fprintf(fp,"\t(diagonal of the mass matrix)\t\t\t(global)\n");
	fprintf(fp,"  Joint X-mass      Y-mass      Z-mass");
	fprintf(fp,"      X-inrta     Y-inrta     Z-inrta\n");
	for (j=1; j <= nJ; j++) {
		k = 6*(j-1);
		fprintf(fp," %5d", j);
		for ( i=1; i<=6; i++ )
			fprintf (fp, " %11.5e", M[k+i][k+i] );
		fprintf(fp,"\n");
	}
	if ( lump )	fprintf(fp,"  Lump masses at joints.\n");
	else		fprintf(fp,"  Use consistent mass matrix.\n");
	fprintf(fp,"N A T U R A L   F R E Q U E N C I E S   & \n");
	fprintf(fp,"M A S S   N O R M A L I Z E D   M O D E   S H A P E S \n");
	fprintf(fp," convergence tolerance: %.3e \n", tol);
	for (m=1; m<=num_modes; m++) { 
	    mpfX = 0.0;	for (i=1; i<=DoF; i++)    mpfX += V[i][m]*msX[i];
	    mpfY = 0.0;	for (i=1; i<=DoF; i++)    mpfY += V[i][m]*msY[i];
	    mpfZ = 0.0;	for (i=1; i<=DoF; i++)    mpfZ += V[i][m]*msZ[i];
	    fprintf(fp,"  MODE %5d:   f= %lf Hz,  T= %lf sec\n",m,f[m],1./f[m]);
	    fprintf(fp,"\t\tX- modal participation factor = %12.4e \n", mpfX);
	    fprintf(fp,"\t\tY- modal participation factor = %12.4e \n", mpfY);
	    fprintf(fp,"\t\tZ- modal participation factor = %12.4e \n", mpfZ);

	    fprintf(fp,"  Joint   X-dsp       Y-dsp       Z-dsp");
	    fprintf(fp,"       X-rot       Y-rot       Z-rot\n");
	    for (j=1; j<= nJ; j++) {
		fprintf(fp," %5d", j);
		for ( i=5; i>=0; i-- )	fprintf (fp, " %11.3e", V[6*j-i][m] );
                fprintf(fp,"\n");
            }
	}

	fprintf(fp,"M A T R I X    I T E R A T I O N S: %d\n", iter );

	fs = sqrt(4.0*PI*PI*f[modes]*f[modes] + tol) / (2.0*PI);

        fprintf(fp,"There are %d modes below %f Hz.", -ok, fs );
        if ( -ok > modes ) {
                fprintf(fp," ... %d modes were not found.\n", -ok-modes );
                fprintf(fp," Try increasing the number of modes in \n");
                fprintf(fp," order to get the missing modes below %f Hz.\n",fs);
        } else  fprintf(fp," ... All %d modes were found.\n", modes );


	free_vector(msX,1,DoF);
	free_vector(msY,1,DoF);
	free_vector(msZ,1,DoF);
	fflush(fp);
	return;
}


/*------------------------------------------------------------------------------
MESH  -  create mesh data of deformed and undeformed mesh, use gnuplot	22feb99
	 useful gnuplot options: set noxtics noytics noztics noborder view nokey
------------------------------------------------------------------------------*/
void mesh(
		char IO_file[], char meshfile[], char plotfile[]
		, char *title, int nJ, int nM, int DoF
		, float *x, float *y, float *z, float *L
		, float *J1, float *J2, float *p, float *D
		, float exg, int anlyz
){
	FILE	*fpmfx, *fpm;
	float	mx, my, mz,	/* coordinates of the member labels	*/
		*vector(); 
	int	j1, j2, i, j, m, X=0, Y=0, Z=0, 
		Strcat(), Strcpy();
	char	meshfl[64], str[10], D3 = '#';
	void	bent_beam(), free_vector();
	void	exit();
        time_t  now;            /* modern time variable type    (DJGPP) */

        (void) time(&now);
 


	Strcpy(meshfl,meshfile);
	str[0]='f'; str[1]='\0';
	Strcat(meshfl,str);

	if ((fpmfx = fopen (meshfl, "w")) == NULL) {
		printf (" error: cannot open meshfile: %s\n", meshfile);
		exit(1);
	}

	if ((fpm = fopen (meshfile, "w")) == NULL) {
		printf (" error: cannot open meshfile: %s\n", meshfile);
		exit(1);
	}

	if (!anlyz) exg = 0.0;



	fprintf(fpm,"# FRAME ANALYSIS RESULTS  http://www.duke.edu/~hpgavin/frame/\n");
	fprintf(fpm,"# %s\n", title );
        fprintf(fpm,"# %s", ctime(&now) );
	fprintf(fpm,"# M E S H   D A T A   (global coordinates)");
	fprintf(fpm," deflection exaggeration: %.1f\n", exg );
	fprintf(fpm,"# Joint      X           Y           Z");
	fprintf(fpm,"          X-dsp       Y-dsp       Z-dsp\n");

	fprintf(fpmfx,"# FRAME ANALYSIS RESULTS  http://www.duke.edu/~hpgavin/frame/\n");
	fprintf(fpmfx,"# %s\n", title );
        fprintf(fpmfx,"# %s", ctime(&now) );
	fprintf(fpmfx,"# F L E X E D   M E S H   D A T A ");
	fprintf(fpmfx,"  deflection exaggeration: %.1f\n", exg );
	fprintf(fpmfx,"#       X-dsp        Y-dsp        Z-dsp\n");

	for (m=1; m<=nM; m++) {

		bent_beam ( fpmfx, J1[m], J2[m], x,y,z, L[m], p[m], D, exg );

		j = J1[m];	i = 6*(j-1);
		fprintf (fpm,"%5d %11.3e %11.3e %11.3e", j, x[j],y[j],z[j]);
		fprintf (fpm," %11.3e %11.3e %11.3e\n",
		x[j]+exg*D[i+1], y[j]+exg*D[i+2], z[j]+exg*D[i+3]);
		j = J2[m];	i = 6*(j-1);
		fprintf (fpm,"%5d %11.3e %11.3e %11.3e", j, x[j],y[j],z[j]);
		fprintf (fpm," %11.3e %11.3e %11.3e\n",
		x[j]+exg*D[i+1], y[j]+exg*D[i+2], z[j]+exg*D[i+3]);
		fprintf(fpm,"\n\n");
	}

	for ( j=1; j<=nJ; j++ ) {
		if (x[j] != 0.0) X=1;	/* check for three-dimensional frame */
		if (y[j] != 0.0) Y=1;
		if (z[j] != 0.0) Z=1;
	}
	if ( X && Y && Z ) D3 = ' ';

	fclose(fpmfx);
	fclose(fpm);

	if ((fpm = fopen (plotfile, "w")) == NULL) {
		printf (" error: cannot open plot file: %s\n", plotfile);
		exit(1);
	}
	fprintf(fpm,"# FRAME ANALYSIS RESULTS  http://www.duke.edu/~hpgavin/frame/\n");
	fprintf(fpm,"# %s\n", title );
	fprintf(fpm,"# %s", ctime(&now) );
	fprintf(fpm,"# M E S H   A N N O T A T I O N   F I L E \n");
	fprintf(fpm,"set title \"%s\\n", title );
	fprintf(fpm,"analysis file: %s ", IO_file );
	fprintf(fpm,"  deflection exaggeration: %.1f\"\n", exg );
	fprintf(fpm,"set autoscale\n");
	fprintf(fpm,"set noborder\n");
	fprintf(fpm,"set pointsize 1.0\n");
	fprintf(fpm,"set xtics; set ytics; set ztics; \n");
	fprintf(fpm,"set nozeroaxis\n");
	fprintf(fpm,"set nokey\n");
	fprintf(fpm,"set nolabel\n");
	
	fprintf(fpm,"# NODE NUMBER LABELS\n");
	for (j=1; j<=nJ; j++) 
		fprintf(fpm,"set label ' %d' at %12.4e, %12.4e, %12.4e\n",
							j, x[j], y[j], z[j] );

	fprintf(fpm,"# MEMBER NUMBER LABELS\n");
	for (m=1; m<=nM; m++) {
		j1 = J1[m];	j2 = J2[m];
		mx = 0.5 * ( x[j1] + x[j2] );
		my = 0.5 * ( y[j1] + y[j2] );
		mz = 0.5 * ( z[j1] + z[j2] );
		fprintf(fpm,"set label ' %d' at %12.4e, %12.4e, %12.4e\n",
								m, mx, my, mz );
	}
	fprintf(fpm,"plot '%s' u 2:3 t 'undeformed mesh' w lp ", meshfile);
	if (!anlyz) fprintf(fpm,"lw 2 lt 1 pt 6 \n");
	else fprintf(fpm,"lw 1 lt 5 pt 6, '%s' u 1:2 t 'deformed mesh' w l lw 2 lt 3\n", meshfl );

	fprintf(fpm,"%c set parametric\n", D3 );
	fprintf(fpm,"%c set view 60, 70, 1 \n", D3 );
	fprintf(fpm,"%c set nokey\n", D3 );
	fprintf(fpm,"%c set xlabel 'x'\n", D3 );
	fprintf(fpm,"%c set ylabel 'y'\n", D3 );
	fprintf(fpm,"%c set zlabel 'z'\n", D3 );
/*	fprintf(fpm,"%c set nolabel\n", D3 );	*/
	fprintf(fpm,"%c splot '%s' u 2:3:4 t 'undeformed mesh' w lp ",
								D3, meshfile );
	if (!anlyz) fprintf(fpm," lw 2 lt 1 pt 6 \n");
	else fprintf(fpm," lw 1 lt 5 pt 6, '%s' u 1:2:3 t 'deformed mesh' w l lw 2 lt 3\n",meshfl);

	fclose(fpm);

	return;
}


/*------------------------------------------------------------------------------
MODAL_MESH  -  create mesh data of the mode-shape meshes, use gnuplot	19oct98
	 useful gnuplot options: set noxtics noytics noztics noborder view nokey
------------------------------------------------------------------------------*/
void modal_mesh(
		char IO_file[], char meshfile[], char modefile[]
		, char plotfile[], char *title
		, int nJ, int nM, int DoF, int modes
		, float *x, float *y, float *z, float *L
		, float *J1, float *J2, float *p
		, float **M, float *f, float **V
		, float exg, int anlyz
){
	FILE	*fpm;
	float	mx, my, mz,	/* coordinates of the member labels	*/
		mpfX, mpfY, mpfZ,	/* mode participation factors	*/
		*msX, *msY, *msZ,
		*v,		/* a mode-shape vector */
		*vector();

	int	i, j, m,n, X=0, Y=0, Z=0,
		Strcat(), Strcpy(); 
	char	D3 = '#', s1[16],  s2[16], modefl[64];
	void	bent_beam(), itoa(), free_vector(), exit();


	msX = vector(1,DoF);
	msY = vector(1,DoF);
	msZ = vector(1,DoF);
	v   = vector(1,DoF);

	for (i=1; i<=DoF; i++) {	/* modal participation factors */
		msX[i] = msY[i] = msZ[i] = 0.0;
		for (j=1; j<=DoF; j+=6) msX[i] += M[i][j];
		for (j=2; j<=DoF; j+=6) msY[i] += M[i][j];
		for (j=3; j<=DoF; j+=6) msZ[i] += M[i][j];
	}
	
	if (!anlyz) exg = 0.0;

	for (m=1; m<=modes; m++) {

	  Strcpy(modefl,modefile);
	  s1[0]='-'; s1[1]='\0'; itoa(m,s2,2);  Strcat(s1,s2);  Strcat(modefl,s1);

	  if ((fpm = fopen (modefl, "w")) == NULL) {
		printf (" error: cannot open modal mesh file: %s\n", modefl);
		exit(1);
	  }

	fprintf(fpm,"# FRAME ANALYSIS RESULTS  http://www.duke.edu/~hpgavin/frame/\n");
	fprintf(fpm,"# %s\n", title );
	  fprintf(fpm,"# M O D E   S H A P E   D A T A   F O R   M O D E");
	  fprintf(fpm,"   %d\t(global coordinates)\n", m );
	  fprintf(fpm,"# deflection exaggeration: %.1f\n\n", exg );
	  mpfX = 0.0;	for (i=1; i<=DoF; i++)    mpfX += V[i][m]*msX[i];
	  mpfY = 0.0;	for (i=1; i<=DoF; i++)    mpfY += V[i][m]*msY[i];
	  mpfZ = 0.0;	for (i=1; i<=DoF; i++)    mpfZ += V[i][m]*msZ[i];
	  fprintf(fpm,"# MODE %5d:   f= %lf Hz, T= %lf sec\n", m,f[m],1./f[m]);
	  fprintf(fpm,"#\t\tX- modal participation factor = %12.4e \n", mpfX);
	  fprintf(fpm,"#\t\tY- modal participation factor = %12.4e \n", mpfY);
	  fprintf(fpm,"#\t\tZ- modal participation factor = %12.4e \n", mpfZ);

	  for (i=1; i<=DoF; i++)	v[i] = V[i][m];

	  fprintf(fpm,"#      X-dsp       Y-dsp       Z-dsp\n\n");

	  for (n=1; n<=nM; n++) 

		bent_beam ( fpm, J1[n], J2[n], x,y,z, L[n], p[n], v, exg );

	  for ( j=1; j<=nJ; j++ ) {
		if (x[j] != 0.0) X=1;	/* check for three-dimensional frame */
		if (y[j] != 0.0) Y=1;
		if (z[j] != 0.0) Z=1;
	  }
	  if ( X && Y && Z ) D3 = ' ';

	  fclose(fpm);

	  if ((fpm = fopen (plotfile, "a")) == NULL) {
		printf (" error: cannot append plot file: %s\n",plotfile);
		exit(1);
	  }
	  fprintf(fpm,"pause -1\n");
	  fprintf(fpm,"set nolabel\n");
	  fprintf(fpm,"set title '%s     mode %d     %lf Hz'\n",IO_file,m,f[m]);
	  fprintf(fpm,"plot '%s' u 2:3 t 'undeformed mesh' w l ", meshfile );
	  if (!anlyz) fprintf(fpm," lw 2 lt 1 \n");
	  else fprintf(fpm," lw 1 lt 5 , '%s' u 1:2 t 'mode-shape %d' w l lw 2 lt 3\n",
								modefl, m );
	  fprintf(fpm,"%c pause -1\n", D3 );
	  fprintf(fpm,"%c set nokey\n", D3 );
	  fprintf(fpm,"%c splot '%s' u 2:3:4 t 'undeformed mesh' w l ",
								D3, meshfile);
	  if (!anlyz) fprintf(fpm," lw 2 lt 1 \n");
	  else fprintf(fpm," lw 1 lt 5 , '%s' u 1:2:3 t 'mode-shape %d' w l lw 2 lt 3\n",
								modefl, m );

	  fclose(fpm);

        }

	free_vector(msX,1,DoF);
	free_vector(msY,1,DoF);
	free_vector(msZ,1,DoF);
	free_vector(v,1,DoF);

	return;
}

/*------------------------------------------------------------------------------
ANIMATE -  create mesh data of animated mode-shape meshes, use gnuplot	16dec98
	 useful gnuplot options: set noxtics noytics noztics noborder view nokey
	 mpeg movie example:   % convert mesh_file-03-f-*.ps mode-03.mpeg
	 ... requires ImageMagick and mpeg2vidcodec packages
------------------------------------------------------------------------------*/
void animate(
	char IO_file[], char meshfile[], char modefile[], char plotfile[]
	, char *title
	, int anim[]
	, int nJ, int nM, int DoF, int modes
	, float *x, float *y, float *z, float *L, float *p
	, int *J1, int *J2, float *f, float **V
	, float exg
	, int pan
){
	FILE	*fpm;
	float	mx, my, mz,	/* coordinates of the member labels	*/
		x_min = 0, x_max = 0,
		y_min = 0, y_max = 0,
		z_min = 0, z_max = 0,
		rot_x_init = 70,	/* inital x-rotation in 3D animation */
		rot_x_final = 60,	/* final  x-rotation in 3D animation */
		rot_z_init = 100,	/* inital z-rotation in 3D animation */
		rot_z_final = 120,	/* final  z-rotation in 3D animation */
		zoom_init = 1.1,	/* inital zoom scale in 3D animation */
		zoom_final = 1.1,	/* final  zoom scale in 3D animation */
		frames = 25,	/* number of frames in animation	*/
		ex=10,		/* an exageration factor, for animation */
		*v, *vector();

	int	fr, i,j, m,n, X=0, Y=0, Z=0, j1,j2, c, CYCLES=3,
		frame_number = 0,
		total_frames,	/* total number of frames in animation */
		Strcat(), Strcpy();
	char	D3 = '#',
		Movie = '#',	/* use '#' for no-movie  -OR-  ' ' for movie */
		s1[16], s2[16], modefl[64], framefl[64];
	void	itoa(), bent_beam(), free_vector(), exit();

	for (j=1; j<=nJ; j++) {		/* check for three-dimensional frame */
		if (x[j] != 0.0) X=1;
		if (y[j] != 0.0) Y=1;
		if (z[j] != 0.0) Z=1;
		if ( x[j] < x_min ) x_min = x[j];
		if ( y[j] < y_min ) y_min = y[j];
		if ( z[j] < z_min ) z_min = z[j];
		if ( x_max < x[j] ) x_max = x[j];
		if ( y_max < y[j] ) y_max = y[j];
		if ( z_max < z[j] ) z_max = z[j];
	}
	if ( X && Y && Z ) D3 = ' ';


	if ((fpm = fopen (plotfile, "a")) == NULL) {
		printf (" error: cannot append plot file: %s\n",plotfile);
		exit(1);
	}
	i = 0;
	while ( (m = anim[i]) != 0 && i < 20) { 
	 if ( i==0 ) {

	   fprintf(fpm,"\n# --- M O D E   S H A P E   A N I M A T I O N ---\n");
	   fprintf(fpm,"set noborder\n");
	   fprintf(fpm,"set nozeroaxis\n");
	   fprintf(fpm,"set autoscale\n");
	   fprintf(fpm,"set noxtics; set noytics; set noztics; \n");
	   fprintf(fpm,"set nokey\n");
	   if ( x_min != x_max )
		fprintf(fpm,"set xrange [ %lf : %lf ] \n",
	 		x_min-0.2*(x_max-x_min), x_max+0.2*(x_max-x_min) );
	   else fprintf(fpm,"set xrange [ %lf : %lf ] \n",
			x_min-exg, x_max+exg );
	   if (y_min != y_max)
		fprintf(fpm,"set yrange [ %lf : %lf ] \n",
	 		y_min-0.2*(y_max-y_min), y_max+0.2*(y_max-y_min) );
	   else fprintf(fpm,"set yrange [ %lf : %lf ] \n",
			y_min-exg, y_max+exg );
	   if (z_min != z_max)
	   	fprintf(fpm,"set zrange [ %lf : %lf ] \n",
			z_min-0.2*(z_max-z_min), z_max+0.2*(z_max-z_min) );
	   else fprintf(fpm,"set zrange [ %lf : %lf ] \n",
			z_min-exg, z_max+exg );
 	
	   fprintf(fpm,"%c set parametric\n", D3 );
	   fprintf(fpm,"%c set view 60, 70, 1 \n", D3 );
	   fprintf(fpm,"%c set xlabel \n", D3 );
	   fprintf(fpm,"%c set ylabel \n", D3 );
	   fprintf(fpm,"%c set zlabel \n", D3 );
	   fprintf(fpm,"%c set nolabel \n", D3 );
	 }

	 fprintf(fpm,"pause -1 \n");
	 fprintf(fpm,"set title '%s     mode %d      %lf Hz'\n",IO_file,m,f[m]);

	 frame_number = 0;
	 total_frames = 2*CYCLES*frames;
	 for ( c=1; c <= CYCLES; c++ ) { 
	  for ( fr=0; fr<=frames; fr++ ) {

	    Strcpy(modefl,modefile);
	    Strcpy(framefl,modefile);
	    s1[0] = '-';  s1[1] = '\0';  itoa(m,s2,2);  Strcat(s1,s2); 
	    Strcat(framefl,s1);
	    Strcat(s1,".");  itoa(fr,s2,3);  Strcat(s1,s2);  Strcat(modefl,s1);
	    s1[0] = '-'; s1[1] = 'f'; s1[2] = '-'; s1[3] = '\0';
	    itoa(frame_number++,s2,3); Strcat(s1,s2); Strcat(framefl,s1);
	    s1[0] = '.'; s1[1] = 'p'; s1[2] = 's'; s1[3] = '\0';
	    Strcat(framefl,s1);

	    if ( D3 == '#' ) {
		fprintf(fpm,"plot '%s' u 2:3 w l lw 1 lt 5, ", meshfile );
	 	fprintf(fpm," '%s' u 1:2 w l lw 2 lt 3 ;", modefl );
	    } else {
	      if (pan) 
 	        fprintf(fpm,"%c set view %5.1f, %5.1f, %4.2f \n", D3,
		rot_x_init + (rot_x_final-rot_x_init)*frame_number/total_frames,
		rot_z_init + (rot_z_final-rot_z_init)*frame_number/total_frames,
		zoom_init + (zoom_final-zoom_init)*frame_number/total_frames );
	      fprintf(fpm,"%c splot '%s' u 2:3:4 w l lw 1 lt 5, ",D3,meshfile);
	      fprintf(fpm," '%s' u 1:2:3 w l lw 2 lt 3;", modefl );
	    } 
	    if ( fr==0 && c==1 )	fprintf(fpm,"  pause 1.5 \n");
	    else			fprintf(fpm,"  pause 0.05 \n");
	    fprintf(fpm,"%c  load 'saveplot';\n",Movie);
	    fprintf(fpm,"%c  !mv my-plot.ps %s\n", Movie, framefl );    
	  }
	  for ( fr = frames-1; fr > 0; fr-- ) {

	    Strcpy(modefl,modefile);
	    Strcpy(framefl,modefile);
	    s1[0] = '-';  s1[1] = '\0';  itoa(m,s2,2);  Strcat(s1,s2); 
	    Strcat(framefl,s1);
	    Strcat(s1,".");  itoa(fr,s2,3);  Strcat(s1,s2);  Strcat(modefl,s1);
	    s1[0] = '-'; s1[1] = 'f'; s1[2] = '-'; s1[3] = '\0';
	    itoa(frame_number++,s2,3); Strcat(s1,s2); Strcat(framefl,s1);
	    s1[0] = '.'; s1[1] = 'p'; s1[2] = 's'; s1[3] = '\0';
	    Strcat(framefl,s1);

	    if ( D3 == '#' ) {
	 	fprintf(fpm,"plot '%s' u 2:3 w l lw 1 lt 5, ", meshfile );
		fprintf(fpm," '%s' u 1:2 w l lw 2 lt 3;", modefl );
	    } else {
	      if (pan)
	        fprintf(fpm,"%c set view %5.1f, %5.1f, %4.2f \n", D3,
		rot_x_init + (rot_x_final-rot_x_init)*frame_number/total_frames,
		rot_z_init + (rot_z_final-rot_z_init)*frame_number/total_frames,
		zoom_init + (zoom_final-zoom_init)*frame_number/total_frames );
	      fprintf(fpm,"%c splot '%s' u 2:3:4 w l lw 1 lt 5, ",D3,meshfile);
	      fprintf(fpm," '%s' u 1:2:3 w l lw 2 lt 3;", modefl );
	    } 
	    fprintf(fpm,"  pause 0.05 \n");
	    fprintf(fpm,"%c  load 'saveplot';\n",Movie);
	    fprintf(fpm,"%c  !mv my-plot.ps %s\n", Movie, framefl );    
	  }
	 }
	 fr = 0;

	 Strcpy(modefl,modefile);
	 s1[0] = '-';  s1[1] = '\0';  itoa(m,s2,2);  Strcat(s1,s2); 
	 Strcat(s1,".");  itoa(fr,s2,3);  Strcat(s1,s2);  Strcat(modefl,s1);

	 if ( D3 == '#' ) {
	 	fprintf(fpm,"plot '%s' u 2:3 w l lw 2 lt 5, ", meshfile );
		fprintf(fpm," '%s' u 1:2 w l lw 3 lt 3 \n", modefl );
	 } else {
		fprintf(fpm,"%c splot '%s' u 2:3:4 w l lw 2 lt 5, ",D3,meshfile);
		fprintf(fpm," '%s' u 1:2:3 w l lw 3 lt 3 \n", modefl );
	 } 

	 i++;
	}
	fclose(fpm);

	v = vector(1,DoF);

	i = 0;
	while ( (m = anim[i]) != 0 ) { 
	  for ( fr=0; fr<=frames; fr++ ) {

	    Strcpy(modefl,modefile);
	    s1[0] = '-';  s1[1] = '\0';  itoa(m,s2,2);  Strcat(s1,s2); 
	    Strcat(s1,".");  itoa(fr,s2,3);  Strcat(s1,s2);  Strcat(modefl,s1);

	    if ((fpm = fopen (modefl, "w")) == NULL) {
		printf (" error: cannot open modal mesh file: %s\n", modefl);
		exit(1);
	    }

	    fprintf(fpm,"# FRAME ANALYSIS RESULTS  http://www.duke.edu/~hpgavin/frame/\n");
	    fprintf(fpm,"# %s\n", title );
	    fprintf(fpm,"# A N I M A T E D   M O D E   S H A P E   D A T A \n");
	    fprintf(fpm,"# deflection exaggeration: %.1f\n", ex );
	    fprintf(fpm,"# MODE %5d: f= %lf Hz  T= %lf sec\n\n",m,f[m],1./f[m]);

	    ex = exg*cos( PI*fr/frames );

	    for (j=1; j<=DoF; j++)	v[j] = V[j][m];

	    fprintf(fpm,"#      X-dsp       Y-dsp       Z-dsp\n\n");

	    for (n=1; n<=nM; n++) 

		bent_beam ( fpm, J1[n], J2[n], x,y,z, L[n], p[n], v, ex );

	    fclose(fpm);
          }
	  i++;
	}
	free_vector(v,1,DoF);
	return;
}


/*------------------------------------------------------------------------------
BENT_BEAM  -  computes cubic deflection functions from beam end deflections
and beam end rotations.  Saves deflected shapes to a file.  These bent shapes
are exact for mode-shapes, and for frames loaded at their joints.	22feb99
------------------------------------------------------------------------------*/
void bent_beam(
		FILE *fp, int j1, int j2
		, float *x, float *y, float *z
		, float L, float p, float *D
		, float exg
){
	float	t1, t2, t3, t4, t5, t6, t7, t8, t9, 	/* coord xfmn	*/
		u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, 
		*a, *b, **A,
		s, v, w, dx, dy, dz,
		*vector(), **matrix();
	int	i1, i2, pd;
	void	lu_dcmp(), free_vector(), free_matrix(), exit();

	A = matrix(1,4,1,4);
	a = vector(1,4);
	b = vector(1,4);

	coord_trans ( x, y, z, L, j1, j2,
				&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p );

	i1 = 6*(j1-1);	i2 = 6*(j2-1);

		/* compute beam end deflections in local coordinates */

	u1  = exg*(t1*D[i1+1] + t2*D[i1+2] + t3*D[i1+3]);
	u2  = exg*(t4*D[i1+1] + t5*D[i1+2] + t6*D[i1+3]);
	u3  = exg*(t7*D[i1+1] + t8*D[i1+2] + t9*D[i1+3]);

	u4  = t1*D[i1+4] + t2*D[i1+5] + t3*D[i1+6];
	u5  = t4*D[i1+4] + t5*D[i1+5] + t6*D[i1+6];
	u6  = t7*D[i1+4] + t8*D[i1+5] + t9*D[i1+6];

	u7  = exg*(t1*D[i2+1] + t2*D[i2+2] + t3*D[i2+3]);
	u8  = exg*(t4*D[i2+1] + t5*D[i2+2] + t6*D[i2+3]);
	u9  = exg*(t7*D[i2+1] + t8*D[i2+2] + t9*D[i2+3]);

	u10 = t1*D[i2+4] + t2*D[i2+5] + t3*D[i2+6];
	u11 = t4*D[i2+4] + t5*D[i2+5] + t6*D[i2+6];
	u12 = t7*D[i2+4] + t8*D[i2+5] + t9*D[i2+6];

		/* curve-fitting problem for a cubic polynomial */

	a[1] = u2;		b[1] = u3;
	a[2] = u8;   		b[2] = u9;
	a[3] = exg*tan(u6);	b[3] = exg*tan(-u5);
	a[4] = exg*tan(u12);	b[4] = exg*tan(-u11);

	u7 += L;
	A[1][1] = 1.0;   A[1][2] = u1;   A[1][3] = u1*u1;   A[1][4] = u1*u1*u1;
	A[2][1] = 1.0;   A[2][2] = u7;   A[2][3] = u7*u7;   A[2][4] = u7*u7*u7;
	A[3][1] = 0.0;   A[3][2] = 1.;   A[3][3] = 2.*u1;   A[3][4] = 3.*u1*u1;
	A[4][1] = 0.0;   A[4][2] = 1.;   A[4][3] = 2.*u7;   A[4][4] = 3.*u7*u7;
	u7 -= L;

	lu_dcmp ( A, 4, a, 1, 1, &pd );		/* solve for cubic coef's */

	if (!pd) {
	 printf(" j1 = %d  j2 = %d  L = %e  u7 = %e \n", j1, j2, L, u7 );
	 exit(1);
	}

	lu_dcmp ( A, 4, b, 0, 1, &pd );		/* solve for cubic coef's */

	for ( s = u1; s <= 1.01*L+u7; s += (L+u7-u1) / 10.0 ) {

			/* deformed shape in local coordinates */
		v = a[1] + a[2]*s + a[3]*s*s + a[4]*s*s*s;
		w = b[1] + b[2]*s + b[3]*s*s + b[4]*s*s*s;
		
			/* deformed shape in global coordinates */
		dx = t1*s + t4*v + t7*w;
		dy = t2*s + t5*v + t8*w;
		dz = t3*s + t6*v + t9*w;

		fprintf (fp," %12.4e %12.4e %12.4e\n",
					x[j1]+dx, y[j1]+dy, z[j1]+dz );
	}
	fprintf(fp,"\n\n");

	free_matrix(A,1,4,1,4);
	free_vector(a,1,4);
	free_vector(b,1,4);

	return;
}


