/*
 FRAME3DD:
 Static and dynamic structural analysis of 2D and 3D frames and trusses with
 elastic and geometric stiffness.
 ---------------------------------------------------------------------------
 http://www.duke.edu/~hpgavin/frame/
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
	Input/output routines for FRAME.

	@note The file format for FRAME is defined in doc/user_manual.html.
*/

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>

#include "common.h"
#include "frame3dd_io.h"
#include "coordtrans.h"
#include "lu_dcmp.h"
#include "nrutil.h"

#ifndef VERSION
# define VERSION "unknown"
#endif

/* #define MASSDATA_DEBUG */

/* forward decls */

static void my_itoa(int n, char s[], int k); /* specialized for portability */

static void getline_no_comment(
	FILE *fp,    /**< pointer to the file from which to read */
	char *s,     /**< pointer to the string to which to write */
	int lim      /**< the longest anticipated line length  */
);

/*------------------------------------------------------------------------------
READ_JOINT_DATA  -  read joint location data				04jan09
------------------------------------------------------------------------------*/
void read_joint_data(
	FILE *fp, int nJ, vec3 *xyz, float *r
){
	int	i, j;

	for (i=1;i<=nJ;i++) {		/* read joint coordinates	*/
		fscanf(fp, "%d", &j );
		if ( j <= 0 || j > nJ ) {
		    fprintf(stderr,"\nERROR: in joint coordinate data, joint number out of range\n");
		    fprintf(stderr,"(joint id is %d <= 0 or > %d)\n", j, nJ);
		    exit(1);
		}
		fscanf(fp, "%lf %lf %lf %f", &xyz[j].x, &xyz[j].y, &xyz[j].z, &r[j]);
		/* fprintf(stderr,"\nj = %d, pos = (%lf, %lf, %lf), r = %f", j, xyz[j].x, xyz[j].y, xyz[j].z, r[j]); */
		r[j] = fabs(r[j]);
	}
	return;
}


/*------------------------------------------------------------------------------
READ_BEAM_DATA  -  read beam property data				04jan09
------------------------------------------------------------------------------*/
void read_beam_data(
	FILE *fp,
	int nJ, int nB, vec3 *xyz, float *r,
	double *L, double *Le,
	int *J1, int *J2,
	float *Ax, float *Asy, float *Asz,
	float *J, float *Iy, float *Iz, float *E, float *G, float *p
){
	int	j1, j2, i, j;

	for (i=1;i<=nB;i++) {		/* read beam properties	*/
		fscanf(fp, "%d", &j );
		if ( j <= 0 || j > nB ) {
		    fprintf(stderr,"  error in beam property data: Beam number out of range  ");
		    fprintf(stderr,"  Beam number: %d  \n", j);
		    exit(1);
		}
		fscanf(fp, "%d %d", &J1[j], &J2[j] );
		if ( J1[j] <= 0 || J1[j] > nJ || J2[j] <= 0 || J2[j] > nJ ) {
		    fprintf(stderr,"  error in beam property data: joint number out of range  ");
		    fprintf(stderr,"  Beam number: %d \n", j);
		    exit(1);
		}
		fscanf(fp, "%f %f %f", &Ax[j], &Asy[j], &Asz[j] );
		fscanf(fp, "%f %f %f", &J[j],  &Iy[j],  &Iz[j] );
		fscanf(fp, "%f %f %f", &E[j], &G[j], &p[j]);

		p[j] = p[j]*PI/180.0;	/* convert from degrees to radians */

		if ( Ax[j] < 0 || Asy[j] < 0 || Asz[j] < 0 ||
		      J[j] < 0 ||  Iy[j] < 0 ||  Iz[j] < 0	) {
		    fprintf(stderr,"  error in beam property data: beam section property < 0  ");
		    fprintf(stderr,"  Beam number: %d  \n", j);
		    exit(1);
		}
		if ( Ax[j] == 0 ) {
		    fprintf(stderr,"  error in beam property data: cross section area is zero   ");
		    fprintf(stderr,"  Beam number: %d  \n", j);
		    exit(1);
		}
		if ( (Asy[j] == 0 || Asz[j] == 0) && G[j] == 0 ) {
		    fprintf(stderr,"  error in beam property data: a shear area and shear modulus are zero   ");
		    fprintf(stderr,"  Beam number: %d  \n", j);
		    exit(1);
		}
		if ( J[j] == 0 ) {
		    fprintf(stderr,"  error in beam property data: torsional moment of inertia is zero   ");
		    fprintf(stderr,"  Beam number: %d  \n", j);
		    exit(1);
		}
		if ( Iy[j] == 0 || Iz[j] == 0 ) {
		    fprintf(stderr,"  error: cross section bending moment of inertia is zero   ");
		    fprintf(stderr,"  Beam number : %d  \n", j);
		    exit(1);
		}
		if ( E[j] <= 0 || G[j] <= 0 ) {
		    fprintf(stderr,"  error : material elastic modulus E or G f is not positive   ");
		    fprintf(stderr,"  Beam number: %d  \n", j);
		    exit(1);
		}
	}
	for (i=1;i<=nB;i++) {		/* calculate beam lengths	*/
		j1 = J1[i];
		j2 = J2[i];

#define SQ(X) ((X)*(X))
		L[i] =	SQ(xyz[j2].x - xyz[j1].x) +
			SQ(xyz[j2].y-xyz[j1].y) +
			SQ(xyz[j2].z-xyz[j1].z);
#undef SQ

		L[i] = sqrt( L[i] );
		Le[i] = L[i] - r[j1] - r[j2];
		if ( j1 == j2 || L[i] == 0.0 ) {
		   fprintf(stderr,
			" Beams must start and stop at different joints\n");
		   fprintf(stderr,
			" beam element %d  J1= %d J2= %d L= %e\n", i, j1,j2, L[i] );
		   fprintf(stderr,
			" Perhaps beam number %d has not been specified. \n", i );
		   exit(1);
		}
		if ( Le[i] <= 0.0 ) {
		   fprintf(stderr, " Joint radii are too large.\n");
		   fprintf(stderr,
			" beam element %d  J1= %d J2= %d L= %e \n", i, j1,j2, L[i] );
		   fprintf(stderr,
			" r1= %e r2= %e Le= %e \n", r[j1], r[j2], Le[i] );
		   exit(1);
		}
	}

	return;
}


/*------------------------------------------------------------------------------
READ_RUN_DATA  -  read information for analysis                   29dec08
------------------------------------------------------------------------------*/
void read_run_data (
	FILE	*fp, filename,
	int	*shear,
	int	*geom,
	char	*meshpath,
	char	*plotpath,
	double	*exagg,
	int	*anlyz
){
	int	i=0, full_len=0, len=0;
	char	mesh_file[96] = "EMPTY_MESH";


	while ( filename[len++] != '\0' ) /* the length of file filename */ ;
        full_len = len;
	while ( filename[len--] != '.' ) /* the last '.' in filename */ ;
	++len;

	for ( i=0; len < full_len; i++,len++ ) ext[i] = tolower(filename[len]);

	fscanf( fp, "%d %d %s %s %lf %d",
			shear, geom, mesh_file, plotpath, exagg, anlyz );

	if (*shear != 0 && *shear != 1) {
	    fprintf(stderr," Rember to specify shear deformations");
	    fprintf(stderr," with a 0 or a 1 after the beam property info.\n");
	    exit(1);
	}

	if (*geom != 0 && *geom != 1) {
	    fprintf(stderr," Rember to specify geometric stiffness");
	    fprintf(stderr," with a 0 or a 1 after the beam property info.\n");
	    exit(1);
	}

	if ( *exagg < 0.0 ) {
	    fprintf(stderr," Remember to specify an exageration");
	    fprintf(stderr," factor greater than zero\n");
	    exit(1);
	}

        output_path(mesh_file,meshpath,FRAME3DD_PATHMAX,NULL);

	return;
}


/*-----------------------------------------------------------------------------
FRAME3DD_GETLINE -  get line into a character string. from K&R         3feb94
-----------------------------------------------------------------------------*/
void frame3dd_getline (
FILE	*fp,
char    *s,
int     lim
){
    int     c=0, i=0;

    while (--lim > 0 && (c=getc(fp)) != EOF && c != '\n' )
            s[i++] = c;
/*      if (c == '\n')  s[i++] = c;	*/
    s[i] = '\0';
    return;
}


/* platform-dependent path sperator character */

#ifdef WIN32
static const char sep = '\\';
#else
static const char sep = '/';
#endif

/*----------------------------------------------------------------------------
TEMP_FILE_LOCATION
return platform-specific temp file locations -- John Pye, Feb 2009
----------------------------------------------------------------------------*/

static const char *temp_dir(){
#ifdef WIN32
	char *tmp;
	tmp = getenv("TEMP");
	if(tmp==NULL){
		fprintf(stderr,"ERROR: %%TEMP%% environment var not found. This "
			"variable needs to be set so that FRAME3DD knows where to put its "
			"temporary files. Set this variable, the re-run FRAME3DD."
                        "The Frame3DD on-line documentation provides help on this.");
		exit(1);
	}
#else
	const char *tmp = "/tmp";
#endif
	return tmp;
}

void temp_file_location(const char *fname, char fullpath[], const int len){
	const char *tmp = temp_dir();
	int res;
	res = snprintf(fullpath,len,"%s%c%s",tmp,sep,fname);
	if(res > len){
		fprintf(stderr,"ERROR: unable to construct temp filename: overflow\n");
		exit(1);
	}
}

/*------------------------------------------------------------------------------
OUTPUT_PATH
return path for output files using either current directory, or FRAME3DD_OUTDIR
if specified. -- John Pye, Feb 2009.
------------------------------------------------------------------------------*/
void output_path(const char *fname, char fullpath[], const int len, const char *default_outdir){
	assert(fname!=NULL);
	int res;
	if(fname[0]==sep){
		/* absolute output path specified */
		res = snprintf(fullpath,len,"%s",fname);
	}else{
		/* fprintf(stderr,"Generating output path for file '%s'\n",fname); */
		const char *outdir;
		outdir = getenv("FRAME3DD_OUTDIR");
		if(outdir==NULL){
			if(default_outdir==NULL){
				outdir = temp_dir();
			}else{
				outdir = default_outdir;
			}
		}
		res = snprintf(fullpath,len,"%s%c%s",outdir,sep,fname);
	}
	if(res > len){
		fprintf(stderr,"ERROR: unable to construct output filename: overflow.\n");
		exit(1);
	}
	/* fprintf(stderr,"Output file path generated: %s\n",fullpath); */
}


/*-----------------------------------------------------------------------------
PARSE_INPUT                                                             7may03
 remove comments from the input file, and write a 'clean' input file
-----------------------------------------------------------------------------*/
void parse_input(FILE *fp, const char *tpath){
	FILE	*fpc;		/* cleaned inout/output file pointer	*/
	char	line[256];

	if ((fpc = fopen (tpath, "w")) == NULL) {
		fprintf (stderr,"ERROR: cannot open file '%s'\n",tpath);
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
	FILE *fp,   /**< pointer to the file from which to read */
	char *s,    /**< pointer to the string to which to write */
	int lim    /**< the longest anticipated line length  */
){
	int     c=0, i=0;

	while (--lim > 0 && (c=getc(fp)) != EOF && c != '\n'
		&& c != '%' && c != '#' && c != ';' && c != '?' ) {
		if (c != ',')
			s[i++] = c;
		else
			s[i++] = ' ';
	/*      if (c == '\n')  s[i++] = c;     */
	}
	s[i] = '\0';
	if (c != '\n')
		while (--lim > 0 && (c=getc(fp)) != EOF && c != '\n')
		/* read the rest of the line, otherwise do nothing */ ;

	if ( c == EOF ) s[0] = EOF;

	return;
}


/*------------------------------------------------------------------------------
READ_REACTION_DATA - Read fixed joint displacement boundary conditions 29dec09
------------------------------------------------------------------------------*/
void read_reaction_data (
	FILE *fp, int DoF, int nJ, int *nR, int *R, int *sumR
){
	int	i,j,l;

	for (i=1; i<=DoF; i++)	R[i] = 0;

	fscanf(fp,"%d", nR );	/* read restrained degrees of freedom */
	printf(" number of joints with reactions ");
	dots(20);
	printf(" nR = %3d",*nR);
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
	    fprintf(stderr,"  %d imposed reactions.\n", *sumR );
	    fprintf(stderr,"  At least 4 reactions are required to support static loads.\n");
	    /*	exit(1); */
	}
	if ( *sumR >= DoF ) {
	    fprintf(stderr,"  error in reaction data:  Fully restrained structure\n");
	    fprintf(stderr,"  %d imposed reactions >= %d degrees of freedom\n",
								*sumR, DoF );
	    exit(1);
	}

	return;
}


/*------------------------------------------------------------------------------
READ_AND_ASSEMBLE_LOADS  -
read load information data, assemble un-restrained load vectors	9sep08
------------------------------------------------------------------------------*/
void read_and_assemble_loads(
		FILE *fp,
		int nJ, int nB, int nL, int DoF,
		vec3 *xyz,
		double *L, double *Le,
		int *J1, int *J2,
		float *Ax, float *Asy, float *Asz,
		float *Iy, float *Iz, float *E, float *G,
		float *p,
		int *R,
		int shear,
		int *nF, int *nW, int *nP, int *nT, int *nD,
		double **Q,
		double **F_mech, double **F_temp,
		float ***W, float ***P, float ***T,
		float **Dp,
		double ***feF_mech, double ***feF_temp
){
	float	hy, hz;			/* section dimensions in local coords */

	double	Nx1, Vy1, Vz1, Mx1, My1, Mz1,	/* fixed end forces */
		Nx2, Vy2, Vz2, Mx2, My2, Mz2,
		Ksy, Ksz, 		/* shear deformatn coefficients	*/
		a, b,				/* point load locations */
		t1, t2, t3, t4, t5, t6, t7, t8, t9;	/* 3D coord Xfrm coef */
	int	i,j,l, lc, n, j1, j2;

	for (j=1; j<=DoF; j++)
		for (lc=1; lc <= nL; lc++)
			F_mech[lc][j] = F_temp[lc][j] = 0.0;
	for (i=1; i<=12; i++)
		for (n=1; n<=nB; n++)
			for (lc=1; lc <= nL; lc++)
				feF_mech[lc][n][i] = feF_temp[lc][n][i] = 0.0;

	for (i=1; i<=DoF; i++)	for (lc=1; lc<=nL; lc++) Dp[lc][i] = 0.0;

	for (i=1;i<=nB;i++)	for(j=1;j<=12;j++)	Q[i][j] = 0.0;

	for (lc = 1; lc <= nL; lc++) {		/* begin load-case loop */

	  printf(" load case %d of %d: \n", lc, nL );

	  fscanf(fp,"%d", &nF[lc] );		/* joint point loads		*/
	  printf("  number of loaded joints ");
	  dots(27);
	  printf(" nF = %3d\n", nF[lc]);
	  for (i=1; i <= nF[lc]; i++) {	/* ! global structural coordinates ! */
		fscanf(fp,"%d", &j);
		if ( j < 1 || j > nJ ) {
		    fprintf(stderr,"  error in joint load data: joint number out of range  ");
		    fprintf(stderr,"  Joint: %d  \n", j);
		    fprintf(stderr,"  Perhaps you did not specify %d joint loads \n", nF[lc] );
		    exit(1);
		}
		for (l=5; l>=0; l--)	fscanf(fp,"%lf", &F_mech[lc][6*j-l] );
		if ( F_mech[lc][6*j-5]==0 && F_mech[lc][6*j-4]==0 && F_mech[lc][6*j-3]==0 && F_mech[lc][6*j-2]==0 && F_mech[lc][6*j-1]==0 && F_mech[lc][6*j]==0 )
		    fprintf(stderr,"   warning: All joint loads applied at joint %d  are zero\n", j );
	  }

	  fscanf(fp,"%d", &nW[lc] );	/* uniform distributed loads	*/
	  printf("  number of uniform distributed loads ");
	  dots(15);
	  printf(" nW = %3d\n", nW[lc]);
	  if ( nW[lc] < 0 || nW[lc] > nB ) {
		fprintf(stderr,"  error: valid ranges for nW is 0 ... %d \n", nB );
		exit(1);
	  }
	  for (i=1; i <= nW[lc]; i++) {	/* ! local element coordinates ! */
		fscanf(fp,"%d", &n );
		if ( n < 1 || n > nB ) {
		    fprintf(stderr,"  error in uniform distributed loads: element number %d is out of range\n",n);
		    exit(1);
		}
		W[lc][i][1] = (double) n;
		for (l=2; l<=4; l++)	fscanf(fp,"%f", &W[lc][i][l] );

		if ( W[lc][i][2]==0 && W[lc][i][3]==0 && W[lc][i][4]==0 )
		    fprintf(stderr,"   warning: All distributed loads applied to beam element %d  are zero\n", (int)W[i][1] );

		Nx1 = Nx2 = W[lc][i][2]*Le[n] / 2.0;
		Vy1 = Vy2 = W[lc][i][3]*Le[n] / 2.0;
		Vz1 = Vz2 = W[lc][i][4]*Le[n] / 2.0;
		Mx1 = Mx2 = 0.0;
		My1 = -W[lc][i][4]*Le[n]*Le[n] / 12.0;	My2 = -My1;
		Mz1 =  W[lc][i][3]*Le[n]*Le[n] / 12.0;	Mz2 = -Mz1;

		/* debugging
		printf("n=%d Vy=%9.2e Vz=%9.2e My=%9.2e Mz=%9.2e\n",
						n, Vy1,Vz1, My1,Mz1 );	*/

		j1 = J1[n];	j2 = J2[n];

		coord_trans ( xyz, L[n], j1, j2,
			&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p[n] );

		/* debugging
		printf("t1=%5.2f t2=%5.2f t3=%5.2f \n", t1, t2, t3 );
                printf("t4=%5.2f t5=%5.2f t6=%5.2f \n", t4, t5, t6 );
                printf("t7=%5.2f t8=%5.2f t9=%5.2f \n", t7, t8, t9 ); */

		/* {F} = [T]'{Q} */
		feF_mech[lc][n][1]  += ( Nx1*t1 + Vy1*t4 + Vz1*t7 );
		feF_mech[lc][n][2]  += ( Nx1*t2 + Vy1*t5 + Vz1*t8 );
		feF_mech[lc][n][3]  += ( Nx1*t3 + Vy1*t6 + Vz1*t9 );
		feF_mech[lc][n][4]  += ( Mx1*t1 + My1*t4 + Mz1*t7 );
		feF_mech[lc][n][5]  += ( Mx1*t2 + My1*t5 + Mz1*t8 );
		feF_mech[lc][n][6]  += ( Mx1*t3 + My1*t6 + Mz1*t9 );

		feF_mech[lc][n][7]  += ( Nx2*t1 + Vy2*t4 + Vz2*t7 );
		feF_mech[lc][n][8]  += ( Nx2*t2 + Vy2*t5 + Vz2*t8 );
		feF_mech[lc][n][9]  += ( Nx2*t3 + Vy2*t6 + Vz2*t9 );
		feF_mech[lc][n][10] += ( Mx2*t1 + My2*t4 + Mz2*t7 );
		feF_mech[lc][n][11] += ( Mx2*t2 + My2*t5 + Mz2*t8 );
		feF_mech[lc][n][12] += ( Mx2*t3 + My2*t6 + Mz2*t9 );

		/* debugging
		printf("n=%d ", n);
		for (l=1;l<=12;l++) {
			if (feF_mech[lc][n][l] != 0)
			   printf(" feF %d = %9.2e ", l, feF_mech[lc][n][l] );
		}
		printf("\n"); */
	  }

	  fscanf(fp,"%d", &nP[lc] );	/* concentrated point loads	*/
	  printf("  number of concentrated beam element point loads ");
	  dots(3);
	  printf(" nP = %3d\n", nP[lc]);
	  if ( nP[lc] < 0 || nP[lc] > nB ) {
		fprintf(stderr,"  error: valid ranges for nP is 0 ... %d \n", nB );
		exit(1);
	  }
	  for (i=1; i <= nP[lc]; i++) {	/* ! local element coordinates ! */
		fscanf(fp,"%d", &n );
		if ( n < 1 || n > nB ) {
		    fprintf(stderr,"   error in internal point loads: beam number %d is out of range\n",n);
		    exit(1);
		}
		P[lc][i][1] = (double) n;
		for (l=2; l<=5; l++)	fscanf(fp,"%f", &P[lc][i][l] );
		a = P[lc][i][5];	b = L[n] - a;

		if ( a < 0 || L[n] < a || b < 0 || L[n] < b ) {
		    fprintf(stderr,"  error in point load data: Point load coord. out of range\n");
		    fprintf(stderr,"  Beam number: %d  L: %lf  load coord.: %lf\n",
							n, L[n], P[lc][i][5] );
		    exit(1);
		}

		if ( shear ) {
			Ksy = G[n]*Asy[n]*Le[n]*Le[n] / (12.*E[n]*Iz[n]);
			Ksz = G[n]*Asz[n]*Le[n]*Le[n] / (12.*E[n]*Iy[n]);
		} else	Ksy = Ksz = 0.0;


		Nx1 = P[lc][i][2]*a/L[n];
		Nx2 = P[lc][i][2]*b/L[n];
		Vy1 = (1./(1.+Ksz))*P[lc][i][3]*b*b*(3.*a + b) / ( L[n]*L[n]*L[n] )+
			(Ksz/(1.+Ksz)) * P[lc][i][3]*b/L[n];
		Vy2 = (1./(1.+Ksz))*P[lc][i][3]*a*a*(3.*b + a) / ( L[n]*L[n]*L[n] )+
			(Ksz/(1.+Ksz)) * P[lc][i][3]*a/L[n];
		Vz1 = (1./(1.+Ksy))*P[lc][i][4]*b*b*(3.*a + b) / ( L[n]*L[n]*L[n] )+
			(Ksy/(1.+Ksy)) * P[lc][i][4]*b/L[n];
		Vz2 = (1./(1.+Ksy))*P[lc][i][4]*a*a*(3.*b + a) / ( L[n]*L[n]*L[n] )+
			(Ksy/(1.+Ksy)) * P[lc][i][4]*a/L[n];
		Mx1 = Mx2 = 0.0;
		My1 = -(1./(1.+Ksy)) * P[lc][i][4]*a*b*b / ( L[n]*L[n] ) -
			(Ksy/(1.+Ksy))* P[lc][i][4]*a*b / (2.*L[n]);
		My2 =  (1./(1.+Ksy)) * P[lc][i][4]*a*a*b / ( L[n]*L[n] ) +
			(Ksy/(1.+Ksy))* P[lc][i][4]*a*b / (2.*L[n]);
		Mz1 =  (1./(1.+Ksz)) * P[lc][i][3]*a*b*b / ( L[n]*L[n] ) +
			(Ksz/(1.+Ksz))* P[lc][i][3]*a*b / (2.*L[n]);
		Mz2 = -(1./(1.+Ksz)) * P[lc][i][3]*a*a*b / ( L[n]*L[n] ) -
			(Ksz/(1.+Ksz))* P[lc][i][3]*a*b / (2.*L[n]);

		j1 = J1[n];	j2 = J2[n];

		coord_trans ( xyz, L[n], j1, j2,
			&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p[n] );

		/* {F} = [T]'{Q} */
		feF_mech[lc][n][1]  += ( Nx1*t1 + Vy1*t4 + Vz1*t7 );
		feF_mech[lc][n][2]  += ( Nx1*t2 + Vy1*t5 + Vz1*t8 );
		feF_mech[lc][n][3]  += ( Nx1*t3 + Vy1*t6 + Vz1*t9 );
		feF_mech[lc][n][4]  += ( Mx1*t1 + My1*t4 + Mz1*t7 );
		feF_mech[lc][n][5]  += ( Mx1*t2 + My1*t5 + Mz1*t8 );
		feF_mech[lc][n][6]  += ( Mx1*t3 + My1*t6 + Mz1*t9 );

		feF_mech[lc][n][7]  += ( Nx2*t1 + Vy2*t4 + Vz2*t7 );
		feF_mech[lc][n][8]  += ( Nx2*t2 + Vy2*t5 + Vz2*t8 );
		feF_mech[lc][n][9]  += ( Nx2*t3 + Vy2*t6 + Vz2*t9 );
		feF_mech[lc][n][10] += ( Mx2*t1 + My2*t4 + Mz2*t7 );
		feF_mech[lc][n][11] += ( Mx2*t2 + My2*t5 + Mz2*t8 );
		feF_mech[lc][n][12] += ( Mx2*t3 + My2*t6 + Mz2*t9 );
	  }

	  fscanf(fp,"%d", &nT[lc] );		/* thermal loads		*/
	  printf("  number of beams with temperature changes ");
	  dots(10);
	  printf(" nT = %3d\n", nT[lc] );
	  if ( nT[lc] < 0 || nT[lc] > nB ) {
		fprintf(stderr,"  error: valid ranges for nT is 0 ... %d \n", nB );
		exit(1);
	  }
	  for (i=1; i <= nT[lc]; i++) {	/* ! element coordinates ! */
		fscanf(fp,"%d", &n );
		if ( n < 1 || n > nB ) {
		    fprintf(stderr,"  error in temperature loads: beam number %d is out of range\n",n);
		    exit(1);
		}
		T[lc][i][1] = (double) n;
		for (l=2; l<=8; l++)	fscanf(fp,"%f", &T[lc][i][l] );
		a  = T[lc][i][2];
		hy = T[lc][i][3];
		hz = T[lc][i][4];

		if ( hy < 0 || hz < 0 ) {
		    fprintf(stderr,"  error in thermal load data: section dimension < 0\n");
		    fprintf(stderr,"  Beam number: %d  hy: %f  hz: %f\n", n,hy,hz);
		    exit(1);
		}

		Nx2 = (a/4.0)*( T[lc][i][5]+T[lc][i][6]+T[lc][i][7]+T[lc][i][8])*E[n]*Ax[n];
		Nx1 = -Nx2;
		Vy1 = Vy2 = Vz1 = Vz2 = 0.0;
		Mx1 = Mx2 = 0.0;
		My1 =  (a/hz)*(T[lc][i][8]-T[lc][i][7])*E[n]*Iy[n];
		My2 = -My1;
		Mz1 =  (a/hy)*(T[lc][i][5]-T[lc][i][6])*E[n]*Iz[n];
		Mz2 = -Mz1;

		j1 = J1[n];	j2 = J2[n];

		coord_trans ( xyz, L[n], j1, j2,
			&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p[n] );

		/* {F} = [T]'{Q} */
		feF_temp[lc][n][1]  += ( Nx1*t1 + Vy1*t4 + Vz1*t7 );
		feF_temp[lc][n][2]  += ( Nx1*t2 + Vy1*t5 + Vz1*t8 );
		feF_temp[lc][n][3]  += ( Nx1*t3 + Vy1*t6 + Vz1*t9 );
		feF_temp[lc][n][4]  += ( Mx1*t1 + My1*t4 + Mz1*t7 );
		feF_temp[lc][n][5]  += ( Mx1*t2 + My1*t5 + Mz1*t8 );
		feF_temp[lc][n][6]  += ( Mx1*t3 + My1*t6 + Mz1*t9 );

		feF_temp[lc][n][7]  += ( Nx2*t1 + Vy2*t4 + Vz2*t7 );
		feF_temp[lc][n][8]  += ( Nx2*t2 + Vy2*t5 + Vz2*t8 );
		feF_temp[lc][n][9]  += ( Nx2*t3 + Vy2*t6 + Vz2*t9 );
		feF_temp[lc][n][10] += ( Mx2*t1 + My2*t4 + Mz2*t7 );
		feF_temp[lc][n][11] += ( Mx2*t2 + My2*t5 + Mz2*t8 );
		feF_temp[lc][n][12] += ( Mx2*t3 + My2*t6 + Mz2*t9 );
	  }

	  for (n=1; n<=nB; n++) {
	     j1 = J1[n];	j2 = J2[n];
	     for (i=1; i<= 6; i++) F_mech[lc][6*j1- 6+i] += feF_mech[lc][n][i];
	     for (i=7; i<=12; i++) F_mech[lc][6*j2-12+i] += feF_mech[lc][n][i];
	     for (i=1; i<= 6; i++) F_temp[lc][6*j1- 6+i] += feF_temp[lc][n][i];
	     for (i=7; i<=12; i++) F_temp[lc][6*j2-12+i] += feF_temp[lc][n][i];
	  }

	  fscanf(fp,"%d", &nD[lc] );	/* read prescribed displacements */
	  printf("  number of joints with prescribed displacements ");
	  dots(4);
	  printf(" nD = %3d\n", nD[lc] );
	  for (i=1; i <= nD[lc]; i++) {
		fscanf(fp,"%d", &j);
		for (l=5; l >=0; l--) {
			fscanf(fp,"%f", &Dp[lc][6*j-l] );
			if ( R[6*j-l] == 0 && Dp[lc][6*j-l] != 0.0 ) {
			    printf(" Initial displacements can be prescribed");
			    printf(" only at restrained coordinates\n");
			    printf(" joint: %d  dof: %d  R: %d\n",
							j, 6-l, R[6*j-l] );
			    exit(1);
			}
		}
	  }

	}					/* end load-case loop */

	return;
}


/*------------------------------------------------------------------------------
READ_MASS_DATA  -  read element densities and extra inertial mass data	16aug01
------------------------------------------------------------------------------*/
void read_mass_data(
		FILE *fp,
		int nJ, int nB, int *nI,
		float *d, float *BMs,
		float *JMs, float *JMx, float *JMy, float *JMz,
		double *L, float *Ax,
		double *total_mass, double *struct_mass,
		int *nM, int *Mmethod, int *lump,
		char modepath[],
		double *tol, double *shift, int anim[], int *pan
){
/*	double	ms = 0.0; */
	int	chk, j, jnt, m, mem, nA;

	char	mode_file[96] = "EMPTY_MODE";

	*total_mass = *struct_mass = 0.0;

	chk = fscanf ( fp, "%d", nM );

	printf(" number of dynamic modes ");
	dots(28);
	printf(" nM = %3d\n", *nM);

	if ( *nM < 1 || chk != 1 ) {
		*nM = 0;
		return;
	}

	fscanf( fp, "%d", Mmethod );

	printf(" modal analysis method ");
	dots(30);
	printf(" %3d ",*Mmethod);
	if ( *Mmethod == 1 ) printf(" (Subspace-Jacobi)\n");
	if ( *Mmethod == 2 ) printf(" (Stodola)\n");


#ifdef MASSDATA_DEBUG
	FILE	*mf;				// mass data file
	mf = fopen("MassData.txt","w");		// open mass data file
	if ((mf = fopen ("MassData.txt", "w")) == NULL) {
	  fprintf (stderr," error: cannot open file 'MassData.txt'\n");
	  exit(1);
	}
	fprintf(mf,"%% structural mass data \n");
	fprintf(mf,"%% element\tAx\t\tlength\t\tdensity\t\tmass \n");
#endif

	fscanf( fp, "%d", lump );
	fscanf( fp, "%s", mode_file );
	fscanf( fp, "%lf", tol );
	fscanf( fp, "%lf", shift );
	for (m=1; m <= nB; m++) {	/* read inertia data	*/
		fscanf(fp, "%d", &mem );
		fscanf(fp, "%f %f", &d[mem], &BMs[mem] );
		*total_mass  += d[mem]*Ax[mem]*L[mem] + BMs[mem];
		*struct_mass += d[mem]*Ax[mem]*L[mem];
#ifdef MASSDATA_DEBUG
		fprintf(mf," %4d\t\t%12.5e\t%12.5e\t%12.5e\t%12.5e \n",
		 mem, Ax[mem], L[mem], d[mem], d[mem]*Ax[mem]*L[mem] );
#endif
	}

#ifdef MASSDATA_DEBUG
	fclose(mf);
#endif

	/* number of joints with extra inertias */
	fscanf(fp,"%d", nI );
	printf(" number of joints with extra lumped inertia ");
        dots(9);
        printf(" nI = %3d\n",*nI);
	for (j=1; j <= *nI; j++) {
		fscanf(fp, "%d", &jnt );
		if ( jnt < 1 || jnt > nJ ) {
	    		fprintf(stderr,"  error in joint load data: joint number out of range  ");
	    		fprintf(stderr,"  Joint: %d  \n", j);
	    		fprintf(stderr,"  Perhaps you did not specify %d extra masses \n", *nI );
	    		exit(1);
		}
		fscanf(fp, "%f %f %f %f",
			&JMs[jnt], &JMx[jnt], &JMy[jnt], &JMz[jnt] );
		*total_mass += JMs[jnt];

		if ( JMs[jnt]==0 && JMx[jnt]==0 && JMy[jnt]==0 && JMz[jnt]==0 )
	    	fprintf(stderr,"  warning: All extra joint inertia at joint %d  are zero\n", jnt );
	}

	for (m=1;m<=nB;m++) {			/* chec inertia data	*/
	    if ( d[m] < 0.0 || BMs[m] < 0.0 || d[m]+BMs[m] <= 0.0 ) {
		fprintf(stderr,"  error: Non-positive mass or density\n");
		fprintf(stderr,"  d[%d]= %f  BMs[%d]= %f\n",m,d[m],m,BMs[m]);
		exit(1);
	    }
	}
/*	for (m=1;m<=nB;m++) ms += BMs[m]; // consistent mass doesn't agree  */
/*	if ( ms > 0.0 )	    *lump = 1;    // with concentrated masses, BMs  */

	printf(" structural mass ");
	dots(36);
	printf("  %12.4e\n",*struct_mass);
	printf(" total mass ");
	dots(41);
	printf("  %12.4e\n",*total_mass);
	fscanf ( fp, "%d", &nA );
	printf(" number of modes to be animated ");
	dots(21);
	printf(" nA = %3d\n",nA);
	if (nA > 20)
	  printf(" nA = %d, only 20 or fewer modes may be animated\n", nA );
	for ( m = 0; m < 20; m++ )	anim[m] = 0;
	for ( m = 0; m < nA; m++ )	fscanf ( fp, "%d", &anim[m] );

	fscanf ( fp, "%d", pan );

	output_path(mode_file,modepath,FRAME3DD_PATHMAX,NULL);


	return;
}


/*------------------------------------------------------------------------------
READ_CONDENSE   -  read matrix condensation information 	        30aug01
------------------------------------------------------------------------------*/
void read_condensation_data (
		FILE *fp,
		int nJ, int nM,
		int *nC, int *Cdof, int *Cmethod, int *q, int *m
){
	int	i,j,k,  chk, **qm;

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
	printf(" nC = %3d\n", *nC );

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
	 if ( (m[i] < 0 || m[i] > nM) && *Cmethod == 3 ) {
	  fprintf(stderr," error in matrix condensation data: \n");
	  fprintf(stderr,"  error: m[%d] = %d \n",i,m[i]);
	  fprintf(stderr,"  The condensed mode number must be between ");
	  fprintf(stderr,"  1 and %d (modes).\n", nM);
	  exit(1);
         }
	}

	free_imatrix(qm,1, *nC, 1,7);
	return;
}


/*------------------------------------------------------------------------------
WRITE_INPUT_DATA  -  save input data					7nov02
------------------------------------------------------------------------------*/
void write_input_data(
	FILE *fp,
	char *title, int nJ, int nB, int nL,
	int *nD, int nR,
	int *nF, int *nW, int *nP, int *nT,
	vec3 *xyz, float *r,
	int *J1, int *J2,
	float *Ax, float *Asy, float *Asz, float *J, float *Iy, float *Iz,
	float *E, float *G, float *p, double **F, float **Dp,
	int *R,
	float ***W, float ***P, float ***T,
	int shear, int anlyz, int geom
){
	int	i,j,n, lc;
        time_t  now;            /* modern time variable type    (DJGPP) */

        (void) time(&now);

	fprintf(fp,"\n");
	for (i=1; i<=80; i++)	fprintf(fp,"_");
  	fprintf(fp,"\nFRAME3DD version: %s ", VERSION );
	fprintf(fp,"              http://frame3dd.sf.net/\n");
	fprintf(fp,"GPL Copyright (C) 1992-2009, Henri P. Gavin \n");
	fprintf(fp,"FRAME3DD is distributed in the hope that it will be useful");
	fprintf(fp," but with no warranty.\n");
	fprintf(fp,"For details see the GNU Public Licence:");
	fprintf(fp," http://www.fsf.org/copyleft/gpl.html\n");
	for (i=1; i<=80; i++)	fprintf(fp,"_"); fprintf(fp,"\n\n");
	fprintf(fp,"%s\n",title);
	fprintf(fp, "%s", ctime(&now) );
	for (i=1; i<=80; i++)	fprintf(fp,"_"); fprintf(fp,"\n");


	fprintf(fp,"%d JOINTS;    %d MEMBERS;    %d LOAD CASES;\n",nJ,nB,nL);
	fprintf(fp,"%d FIXED JOINTS;   ", nR);

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
		i, xyz[i].x, xyz[i].y, xyz[i].z, r[i],
			R[j+1], R[j+2], R[j+3], R[j+4], R[j+5], R[j+6] );
	}
	fprintf(fp,"B E A M   D A T A\t\t\t\t\t\t\t(local)\n");
	fprintf(fp,"  Beam   J1    J2     Ax   Asy   Asz    ");
	fprintf(fp,"Jxx     Iyy     Izz       E       G roll\n");
	for (i=1; i<= nB; i++) {
		fprintf(fp,"%5d %5d %5d %6.1f %5.1f %5.1f",
					i, J1[i],J2[i], Ax[i], Asy[i], Asz[i] );
		fprintf(fp," %6.1f %7.1f %7.1f %8.1f %7.1f %3.0f\n",
				J[i], Iy[i], Iz[i], E[i], G[i], p[i]*180.0/PI );
	}
	if ( shear )	fprintf(fp,"  Include shear deformations.\n");
	else		fprintf(fp,"  Neglect shear deformations.\n");
	if ( geom )	fprintf(fp,"  Include geometric stiffness.\n");
	else		fprintf(fp,"  Neglect geometric stiffness.\n");

	for (lc = 1; lc <= nL; lc++) {		/* start load case loop */

	  fprintf(fp,"\nL O A D   C A S E   %d   O F   %d  ... \n\n", lc,nL);
	  fprintf(fp," %3d joints  with concentrated loads\n", nF[lc] );
	  fprintf(fp," %3d elements with uniformly distributed loads\n", nW[lc]);
	  fprintf(fp," %3d elements with concentrated point loads\n", nP[lc] );
	  fprintf(fp," %3d elements with temperature loads\n", nT[lc] );
	  fprintf(fp," %3d joints  with prescribed displacements\n", nD[lc] );
	  if ( nF[lc] > 0 || nW[lc] > 0 || nP[lc] > 0 || nT[lc] > 0 ) {
	    fprintf(fp," J O I N T   L O A D S");
	    fprintf(fp,"  +  E Q U I V A L E N T   J O I N T   L O A D S  (global)\n");
	    fprintf(fp,"  Joint       Fx          Fy          Fz");
	    fprintf(fp,"          Mxx         Myy         Mzz\n");
	    for (j=1; j<=nJ; j++) {
		i = 6*(j-1);
		if ( F[lc][i+1]!=0.0 || F[lc][i+2]!=0.0 || F[lc][i+3]!=0.0 ||
		     F[lc][i+4]!=0.0 || F[lc][i+5]!=0.0 || F[lc][i+6]!=0.0 ) {
			fprintf(fp, " %5d", j);
			for (i=5; i>=0; i--) fprintf(fp, " %11.3f", F[lc][6*j-i] );
			fprintf(fp, "\n");
		}
	    }
	  }

	  if ( nW[lc] > 0 ) {
	    fprintf(fp," U N I F O R M   B E A M   L O A D S");
	    fprintf(fp,"\t\t\t\t\t(local)\n");
	    fprintf(fp,"  Beam        Wx               Wy               Wz\n");
	    for (n=1; n<=nW[lc]; n++) {
		fprintf(fp, " %5d", (int) (W[lc][n][1]) );
		for (i=2; i<=4; i++) fprintf(fp, " %16.8f", W[lc][n][i] );
		fprintf(fp, "\n");
	    }
	  }

	  if ( nP[lc] > 0 ) {
	    fprintf(fp," C O N C E T R A T E D   P O I N T   L O A D S");
	    fprintf(fp,"\t\t\t\t(local)\n");
	    fprintf(fp,"  Beam        Px          Py          Pz          x\n");
	    for (n=1; n<=nP[lc]; n++) {
		fprintf(fp, " %5d", (int) (P[lc][n][1]) );
		for (i=2; i<=5; i++) fprintf(fp, " %11.3f", P[lc][n][i] );
		fprintf(fp, "\n");
	    }
	  }

	  if ( nT[lc] > 0 ) {
	    fprintf(fp," B E A M   T E M P E R A T U R E   C H A N G E S");
	    fprintf(fp,"\t\t\t(local)\n");
	    fprintf(fp,"  Beam      coef      hy        hz");
	    fprintf(fp,"        Ty+       Ty-       Tz+       Tz-\n");
	    for (n=1; n<=nT[lc]; n++) {
		fprintf(fp, " %5d", (int) (T[lc][n][1]) );
		fprintf(fp, " %9.2e", T[lc][n][2] );
		for (i=3; i<=8; i++) fprintf(fp, " %9.3f", T[lc][n][i] );
		fprintf(fp, "\n");
	    }
	  }

	  if ( nD[lc] > 0 ) {
	    fprintf(fp,"\nP R E S C R I B E D   D I S P L A C E M E N T S");
	    fprintf(fp,"                        (global)\n");
	    fprintf(fp,"  Joint       Dx          Dy          Dz");
	    fprintf(fp,"          Dxx         Dyy         Dzz\n");
	    for (j=1; j<=nJ; j++) {
		i = 6*(j-1);
		if ( Dp[lc][i+1]!=0.0 || Dp[lc][i+2]!=0.0 || Dp[lc][i+3]!=0.0 ||
		     Dp[lc][i+4]!=0.0 || Dp[lc][i+5]!=0.0 || Dp[lc][i+6]!=0.0 ){
			fprintf(fp, " %5d", j);
			for (i=5; i>=0; i--) fprintf(fp, " %11.3f",
							Dp[lc][6*j-i] );
			fprintf(fp, "\n");
		}
	    }
	  }

	}					/* end load case loop	*/

	if (anlyz) {
	 fprintf(fp,"\nE L A S T I C   S T I F F N E S S   A N A L Y S I S");
	 fprintf(fp,"   via  L D L'  decomposition\n\n");
	}
	else		fprintf(fp,"D A T A   C H E C K   O N L Y\n");
	fflush(fp);
	return;
}


/*------------------------------------------------------------------------------
WRITE_STATIC_RESULTS -  save joint displacements and beam end forces	9sep08
------------------------------------------------------------------------------*/
void write_static_results (
		FILE *fp,
		int nJ, int nB, int nL, int lc, int DoF,
		int *J1, int *J2,
		double *F, double *D, int *R, double **Q,
		double err, int ok
){
	double	disp;
	int	i,j,n;

	if ( ok < 0 ) {
	 fprintf(fp,"  * The Stiffness Matrix is not positive-definite *\n");
	 fprintf(fp,"    Check that all six rigid-body translations are restrained\n");
	 fprintf(fp,"    If geometric stiffness is included, reduce the loads.\n");
/*	 return; */
	}

        fprintf(fp,"\nL O A D   C A S E   %d   O F   %d  ... \n\n", lc, nL);

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
                                fprintf (fp, "    0.0     ");
                        else    fprintf (fp, " %11.6f",  D[6*j-i] );
		}
		fprintf(fp,"\n");
	    }
	}
	fprintf(fp,"B E A M   E N D   F O R C E S");
	fprintf(fp,"\t\t\t\t\t\t(local)\n");
	fprintf(fp,"  Beam   Joint      Nx          Vy         Vz");
	fprintf(fp,"         Txx        Myy        Mzz\n");
	for (n=1; n<= nB; n++) {
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
WRITE_STATIC_CSV -  save joint displacements and beam end forces	31dec08
------------------------------------------------------------------------------*/
void write_static_csv(
		char *argv[],
		char *title,
		int nJ, int nB, int nL, int lc, int DoF,
		int *J1, int *J2,
		double *F, double *D, int *R, double **Q,
		double err, int ok
){

	FILE	*fpcsv;
	int	i,j,n;
	char	*wa;
	char	IOfilename[128];
        time_t  now;            /* modern time variable type    (DJGPP) */

        (void) time(&now);

	i=0;
	j=0;
	while (i<128) {
		IOfilename[j] = argv[1][i];
		if ( IOfilename[j] == '+' ||
		     IOfilename[j] == '-' ||
		     IOfilename[j] == '*' ||
		     IOfilename[j] == '^' ||
                     IOfilename[j] == '.' ||
                     IOfilename[j] == '\0') {
			IOfilename[j] = '_';
			break;
		}
		i++;
		j++;
	}
	IOfilename[++j] = '\0';
	strcat(IOfilename,"out.CSV");


	wa  = "a";
	if (lc == 1) wa = "w";

	if ((fpcsv = fopen (IOfilename, wa)) == NULL) {
	  fprintf (stderr," error: cannot open file %s\n", IOfilename );
	  exit(1);
	}


	if ( lc == 1 ) {
  	 fprintf(fpcsv,"\" FRAME3DD version: %s ", VERSION );
	 fprintf(fpcsv,"              http://frame3dd.sf.net/\"\n");
	 fprintf(fpcsv,"\"GPL Copyright (C) 1992-2009, Henri P. Gavin \"\n");
	 fprintf(fpcsv,"\"FRAME3DD is distributed in the hope that it will be useful");
	 fprintf(fpcsv," but with no warranty.\"\n");
	 fprintf(fpcsv,"\"For details see the GNU Public Licence:");
	 fprintf(fpcsv," http://www.fsf.org/copyleft/gpl.html\"\n");
	 fprintf(fpcsv,"\" %s \"\n",title);
	 fprintf(fpcsv,"\" %s \"\n", ctime(&now) );

	 fprintf(fpcsv,"\" .CSV formatted results of Frame3DD analysis \"\n");
	 fprintf(fpcsv,"\n , Load Case , Displacements , End Forces , Reactions \n");
	 for (i = 1; i <= nL; i++) {
	 	fprintf(fpcsv," First Row , %d , %d , %d , %d  \n",
			i,
			15+(i-1)*(nJ*2+nB*2+10) + 2*nL,
			17+(i-1)*(nJ*2+nB*2+10) + 2*nL+ nJ,
			19+(i-1)*(nJ*2+nB*2+10) + 2*nL+ nJ + 2*nB );
	 	fprintf(fpcsv," Last Row , %d , %d , %d , %d  \n",
			i,
			15+(i-1)*(nJ*2+nB*2+10) + 2*nL + nJ - 1,
			17+(i-1)*(nJ*2+nB*2+10) + 2*nL + nJ + 2*nB - 1,
			19+(i-1)*(nJ*2+nB*2+10) + 2*nL + 2*nJ + 2*nB - 1 );
	 }

	}


	if ( ok < 0 ) {
	 fprintf(fpcsv,"\"  * The Stiffness Matrix is not positive-definite * \"\n");
	 fprintf(fpcsv,"\" Check that all six rigid-body translations are restrained\"\n");
	 fprintf(fpcsv,"\" If geometric stiffness is included, reduce the loads.\"\n");
/*	 return; */
	}


        fprintf(fpcsv,"\n\"L O A D   C A S E   %d   O F   %d  ... \"\n\n", lc, nL);

	fprintf(fpcsv,"\"J O I N T   D I S P L A C E M E N T S");
	fprintf(fpcsv,"  (global)\"\n");
	fprintf(fpcsv,"Joint ,  X-dsp   ,   Y-dsp  ,    Z-dsp");
	fprintf(fpcsv," ,     X-rot  ,    Y-rot   ,   Z-rot\n");
	for (j=1; j<= nJ; j++) {
		fprintf(fpcsv," %5d,", j);
		for ( i=5; i>=0; i-- ) {
                        if ( fabs(D[6*j-i]) < 1.e-8 )
                                fprintf (fpcsv, "    0.0,    ");
                        else    fprintf (fpcsv, " %12.5e,",  D[6*j-i] );
		}
		fprintf(fpcsv,"\n");
	}
	fprintf(fpcsv,"\"B E A M   E N D   F O R C E S");
	fprintf(fpcsv,"  (local)\"\n");
	fprintf(fpcsv,"Beam  , Joint ,    Nx     ,    Vy   ,     Vz");
	fprintf(fpcsv,"   ,     Txx   ,    Myy  ,     Mzz\n");
	for (n=1; n<= nB; n++) {
		fprintf(fpcsv," %5d, %5d,", n, J1[n]);
		if ( fabs(Q[n][1]) < 0.0001 )
			fprintf (fpcsv, "      0.0,  ");
		else    fprintf (fpcsv, " %12.5e,", Q[n][1] );
		for (i=2; i<=6; i++) {
			if ( fabs(Q[n][i]) < 0.0001 )
				fprintf (fpcsv, "      0.0, ");
			else    fprintf (fpcsv, " %12.5e,", Q[n][i] );
                }
		fprintf(fpcsv,"\n");
		fprintf(fpcsv," %5d, %5d,", n, J2[n]);
		if ( fabs(Q[n][7]) < 0.0001 )
			fprintf (fpcsv, "      0.0,  ");
		else    fprintf (fpcsv, " %12.5e,", Q[n][7] );
		for (i=8; i<=12; i++) {
			if ( fabs(Q[n][i]) < 0.0001 )
				fprintf (fpcsv, "      0.0, ");
			else    fprintf (fpcsv, " %12.5e,", Q[n][i] );
		}
		fprintf(fpcsv,"\n");
	}
	fprintf(fpcsv,"\"R E A C T I O N S  (global)\"\n");
	fprintf(fpcsv," Joint  ,    Fx      ,   Fy   ,      Fz");
	fprintf(fpcsv,"   ,     Mxx    ,    Myy    ,    Mzz\n");
	for (j=1; j<=nJ; j++) {
		i = 6*(j-1);
		fprintf(fpcsv, " %5d,", j);
                for (i=5; i>=0; i--) {
                	if ( !R[6*j-i] || fabs(F[6*j-i]) < 0.0001 )
                        	fprintf (fpcsv, "       0.0, ");
			else    fprintf (fpcsv, " %12.5e,", -F[6*j-i] );
		}
		fprintf(fpcsv, "\n");
	}
	fprintf(fpcsv,"\"R M S   E Q U I L I B R I U M    E R R O R:\", %9.3e\n", err );

	fclose(fpcsv);

	return;
}

/*------------------------------------------------------------------------------
WRITE_STATIC_MFILE -  	9sep08
save joint displacements and beam end forces in an m-file
this function interacts with frame_3dd.m, an m-file interface to frame3dd
------------------------------------------------------------------------------*/
void write_static_mfile (
		char *argv[], char *title,
		int nJ, int nB, int nL, int lc, int DoF,
		int *J1, int *J2,
		double *F, double *D, int *R, double **Q,
		double err, int ok
){
	FILE	*fpm;
	int	i,j,n;
	char	*wa;
	char	IOfilename[128];
        time_t  now;            /* modern time variable type    (DJGPP) */

        (void) time(&now);

	i=0;
	j=0;
	while (i<128) {
		IOfilename[j] = argv[1][i];
		if ( IOfilename[j] == '+' ||
		     IOfilename[j] == '-' ||
		     IOfilename[j] == '*' ||
		     IOfilename[j] == '^' ||
                     IOfilename[j] == '.' ||
                     IOfilename[j] == '\0') {
			IOfilename[j] = '_';
			break;
		}
		i++;
		j++;
	}
	IOfilename[++j] = '\0';
	strcat(IOfilename,"out.m");

	wa  = "a";
	if (lc == 1) wa = "w";

	if ((fpm = fopen (IOfilename, wa)) == NULL) {
	  fprintf (stderr," error: cannot open file %s\n", IOfilename );
	  exit(1);
	}

	if ( lc == 1 ) {
  	 fprintf(fpm,"%% FRAME3DD version: %s ", VERSION );
	 fprintf(fpm,"              http://frame3dd.sf.net/\n");
	 fprintf(fpm,"%%GPL Copyright (C) 1992-2009, Henri P. Gavin \n");
	 fprintf(fpm,"%%FRAME3DD is distributed in the hope that it will be useful");
	 fprintf(fpm," but with no warranty.\n");
	 fprintf(fpm,"%%For details see the GNU Public Licence:");
	 fprintf(fpm," http://www.fsf.org/copyleft/gpl.html\n");
	 fprintf(fpm,"%% %s\n",title);
	 fprintf(fpm, "%% %s", ctime(&now) );

	 fprintf(fpm,"%% m-file formatted results of frame3dd analysis\n");
	 fprintf(fpm,"%% to be read by frame_3dd.m\n");
	}


	if ( ok < 0 ) {
	 fprintf(fpm,"%%  The Stiffness Matrix is not positive-definite *\n");
	 fprintf(fpm,"%%  Check that all six rigid-body translations are restrained\n");
	 fprintf(fpm,"%%  If geometric stiffness is included, reduce the loads.\n");
/*	 return; */
	}

        fprintf(fpm,"\n%% L O A D   C A S E   %d   O F   %d  ... \n\n", lc, nL);

	fprintf(fpm,"%% J O I N T   D I S P L A C E M E N T S");
	fprintf(fpm,"\t\t(global)\n");
	fprintf(fpm,"%%\tX-dsp\t\tY-dsp\t\tZ-dsp\t\tX-rot\t\tY-rot\t\tZ-rot\n");
        fprintf(fpm,"D%d=[",lc);
	for (j=1; j<= nJ; j++) {
		for ( i=5; i>=0; i-- ) {
                        if ( fabs(D[6*j-i]) < 1.e-8 )
                                fprintf (fpm, "\t0.0\t");
                        else    fprintf (fpm, "\t%13.6e",  D[6*j-i] );
		}
		if ( j < nJ )	fprintf(fpm," ; \n");
		else		fprintf(fpm," ]'; \n\n");
	}

	fprintf(fpm,"%% B E A M   E N D   F O R C E S");
	fprintf(fpm,"\t\t\t\t(local)\n");
	fprintf(fpm,"%%\tNx_1\t\tVy_1\t\tVz_1\t\tTxx_1\t\tMyy_1\t\tMzz_1\t");
	fprintf(fpm,"  \tNx_2\t\tVy_2\t\tVz_2\t\tTxx_2\t\tMyy_2\t\tMzz_2\n");
        fprintf(fpm,"F%d=[",lc);
	for (n=1; n<= nB; n++) {
		if ( fabs(Q[n][1]) < 0.0001 )
			fprintf (fpm, "\t0.0\t");
		else    fprintf (fpm, "\t%13.6e", Q[n][1] );
		for (i=2; i<=6; i++) {
			if ( fabs(Q[n][i]) < 0.0001 )
				fprintf (fpm, "\t0.0\t");
			else    fprintf (fpm, "\t%13.6e", Q[n][i] );
		}
		if ( fabs(Q[n][7]) < 0.0001 )
			fprintf (fpm, "\t0.0\t");
		else    fprintf (fpm, "\t%13.6e", Q[n][7] );
		for (i=8; i<=12; i++) {
			if ( fabs(Q[n][i]) < 0.0001 )
				fprintf (fpm, "\t0.0\t");
			else    fprintf (fpm, "\t%13.6e", Q[n][i] );
		}
		if ( n < nB )	fprintf(fpm," ; \n");
		else		fprintf(fpm," ]'; \n\n");
	}

	fprintf(fpm,"%% R E A C T I O N S\t\t\t\t(global)\n");
	fprintf(fpm,"%%\tFx\t\tFy\t\tFz\t\tMxx\t\tMyy\t\tMzz\n");
        fprintf(fpm,"R%d=[",lc);
	for (j=1; j<=nJ; j++) {
		i = 6*(j-1);
		for (i=5; i>=0; i--) {
			if ( !R[6*j-i] || fabs(F[6*j-i]) < 0.0001 )
				fprintf (fpm, "\t0.0\t");
			else    fprintf (fpm, "\t%13.6e", -F[6*j-i] );
		}
		if ( j < nJ )	fprintf(fpm," ; \n");
		else		fprintf(fpm," ]'; \n\n");
	}

	fprintf(fpm,"%% R M S   E Q U I L I B R I U M    E R R O R: %9.3e\n", err );
	fprintf(fpm,"\n\n  load Ks \n\n");

	fclose(fpm);

	return;
}


/*------------------------------------------------------------------------------
WRITE_MODAL_RESULTS -  save modal frequencies and mode shapes		16aug01
------------------------------------------------------------------------------*/
void write_modal_results(
		FILE *fp,
		int nJ, int nB, int nI, int DoF,
		double **M, double *f, double **V,
		double total_mass, double struct_mass,
		int iter, int sumR, int nM,
		double shift, int lump, double tol, int ok
){
	int	i, j, k, m, num_modes;
	double	mpfX, mpfY, mpfZ,	/* mode participation factors	*/
		*msX, *msY, *msZ;
	double	fs;

	msX = dvector(1,DoF);
	msY = dvector(1,DoF);
	msZ = dvector(1,DoF);

	for (i=1; i<=DoF; i++) {
		msX[i] = msY[i] = msZ[i] = 0.0;
		for (j=1; j<=DoF; j+=6) msX[i] += M[i][j];
		for (j=2; j<=DoF; j+=6) msY[i] += M[i][j];
		for (j=3; j<=DoF; j+=6) msZ[i] += M[i][j];
	}

	if ( (DoF - sumR) > nM )	num_modes = nM;
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

	fs = sqrt(4.0*PI*PI*f[nM]*f[nM] + tol) / (2.0*PI);

        fprintf(fp,"There are %d modes below %f Hz.", -ok, fs );
        if ( -ok > nM ) {
                fprintf(fp," ... %d modes were not found.\n", -ok-nM );
                fprintf(fp," Try increasing the number of modes in \n");
                fprintf(fp," order to get the missing modes below %f Hz.\n",fs);
        } else  fprintf(fp," ... All %d modes were found.\n", nM );


	free_dvector(msX,1,DoF);
	free_dvector(msY,1,DoF);
	free_dvector(msZ,1,DoF);
	fflush(fp);
	return;
}


/*------------------------------------------------------------------------------
STATIC_MESH  -
create mesh data of deformed and undeformed mesh, use gnuplot	22feb99
useful gnuplot options: set noxtics noytics noztics noborder view nokey
------------------------------------------------------------------------------*/
void static_mesh(
		char IO_file[], char meshpath[], char plotpath[],
		char *title, int nJ, int nB, int nL, int lc, int DoF,
		vec3 *xyz, double *L,
		int *J1, int *J2, float *p, double *D,
		double exagg, int anlyz
){
	FILE	*fpmfx, *fpm;
	double	mx, my, mz;	/* coordinates of the beam number labels */
	int	j1, j2, i, j, m, X=0, Y=0, Z=0;
	char	meshfl[64], str[8], D3 = '#';
	time_t  now;            /* modern time variable type    (DJGPP) */

	strcpy(meshfl,meshpath);
	str[0]='f'; str[1]='\0';	strcat(meshfl,str);
	str[0]='.'; str[1]='\0';	strcat(meshfl,str);
	my_itoa(lc,str,3);		strcat(meshfl,str);

	if ((fpmfx = fopen (meshfl, "w")) == NULL) {
		printf (" error: cannot open meshpath: %s\n", meshpath);
		exit(1);
	}

	if ((fpm = fopen (meshpath, "w")) == NULL) {
		printf (" error: cannot open meshpath: %s\n", meshpath);
		exit(1);
	}

	if (!anlyz) exagg = 0.0;



	fprintf(fpm,"# FRAME3DD ANALYSIS RESULTS  http://frame3dd.sf.net/\n");
	fprintf(fpm,"# %s\n", title );
	fprintf(fpm,"# L O A D  C A S E   %d  of   %d \n", lc, nL );
        fprintf(fpm,"# %s", ctime(&now) );
	fprintf(fpm,"# M E S H   D A T A   (global coordinates)");
	fprintf(fpm," deflection exaggeration: %.1f\n", exagg );
	fprintf(fpm,"# Joint      X           Y           Z");
	fprintf(fpm,"          X-dsp       Y-dsp       Z-dsp\n");

	fprintf(fpmfx,"# FRAME3DD ANALYSIS RESULTS  http://frame3dd.sf.net/\n");
	fprintf(fpmfx,"# %s\n", title );
        fprintf(fpmfx,"# %s", ctime(&now) );
	fprintf(fpmfx,"# F L E X E D   M E S H   D A T A ");
	fprintf(fpmfx,"  deflection exaggeration: %.1f\n", exagg );
	fprintf(fpmfx,"#       X-dsp        Y-dsp        Z-dsp\n");

	for (m=1; m<=nB; m++) {

		bent_beam ( fpmfx, J1[m], J2[m], xyz, L[m], p[m], D, exagg );

		j = J1[m];	i = 6*(j-1);
		fprintf (fpm,"%5d %11.3e %11.3e %11.3e", j, xyz[j].x,xyz[j].y,xyz[j].z);
		fprintf (fpm," %11.3e %11.3e %11.3e\n",
			xyz[j].x + exagg*D[i+1],
			xyz[j].y + exagg*D[i+2],
			xyz[j].z + exagg*D[i+3]
		);

		j = J2[m];	i = 6*(j-1);
		fprintf (fpm,"%5d %11.3e %11.3e %11.3e", j, xyz[j].x,xyz[j].y,xyz[j].z);
		fprintf (fpm," %11.3e %11.3e %11.3e\n",
			xyz[j].x + exagg*D[i+1],
			xyz[j].y + exagg*D[i+2],
			xyz[j].z + exagg*D[i+3]
		);
		fprintf(fpm,"\n\n");
	}

	for ( j=1; j<=nJ; j++ ) {
		if (xyz[j].x != 0.0) X=1;	/* check for three-dimensional frame */
		if (xyz[j].y != 0.0) Y=1;
		if (xyz[j].z != 0.0) Z=1;
	}
	if ( X && Y && Z ) D3 = ' ';

	fclose(fpmfx);
	fclose(fpm);

	if (lc == 1) {
	    if ((fpm = fopen (plotpath, "w")) == NULL) {
		printf (" error: cannot open plot file: %s\n", plotpath);
		exit(1);
	    }
	} else {
	    if ((fpm = fopen (plotpath, "a")) == NULL) {
		printf (" error: cannot open plot file: %s\n", plotpath);
		exit(1);
	    }
	}

	if (lc == 1) {		/* first load case */

	 fprintf(fpm,"# FRAME3DD ANALYSIS RESULTS  http://frame3dd.sf.net/\n");
	 fprintf(fpm,"# %s\n", title );
	 fprintf(fpm,"# %s", ctime(&now) );
	 fprintf(fpm,"# M E S H   A N N O T A T I O N   F I L E \n");

	 fprintf(fpm,"set title \"%s\\n", title );
	 fprintf(fpm,"analysis file: %s ", IO_file );
	 fprintf(fpm,"  deflection exaggeration: %.1f ", exagg );
	 fprintf(fpm,"  load case %d of %d \"\n", lc, nL );

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
					j, xyz[j].x,xyz[j].y,xyz[j].z );

	 fprintf(fpm,"# MEMBER NUMBER LABELS\n");
	 for (m=1; m<=nB; m++) {
		j1 = J1[m];	j2 = J2[m];
		mx = 0.5 * ( xyz[j1].x + xyz[j2].x );
		my = 0.5 * ( xyz[j1].y + xyz[j2].y );
		mz = 0.5 * ( xyz[j1].z + xyz[j2].z );
		fprintf(fpm,"set label ' %d' at %12.4e, %12.4e, %12.4e\n",
								m, mx, my, mz );
	 }
	 fprintf(fpm,"plot '%s' u 2:3 t 'undeformed mesh' w lp ", meshpath);
	 if (!anlyz) fprintf(fpm,"lw 2 lt 1 pt 6 \n");
	 else fprintf(fpm,"lw 1 lt 5 pt 6, '%s' u 1:2 t 'load case %d of %d' w l lw 2 lt 3\n", meshfl, lc, nL );

	 fprintf(fpm,"%c set parametric\n", D3 );
	 fprintf(fpm,"%c set view 60, 70, 1 \n", D3 );
	 fprintf(fpm,"%c set nokey\n", D3 );
	 fprintf(fpm,"%c set xlabel 'x'\n", D3 );
	 fprintf(fpm,"%c set ylabel 'y'\n", D3 );
	 fprintf(fpm,"%c set zlabel 'z'\n", D3 );
/*	 fprintf(fpm,"%c set nolabel\n", D3 );	*/
	 fprintf(fpm,"%c splot '%s' u 2:3:4 t 'load case %d of %d' w lp ",
							D3, meshpath, lc, nL );
	 if (!anlyz) fprintf(fpm," lw 2 lt 1 pt 6 \n");
	 else fprintf(fpm," lw 1 lt 5 pt 6, '%s' u 1:2:3 t 'load case %d of %d' w l lw 2 lt 3\n",meshfl, lc, nL );

	} else { 		/* additional load cases */

	 fprintf(fpm,"pause -1\n");

	 fprintf(fpm,"set title \"%s\\n", title );
	 fprintf(fpm,"analysis file: %s ", IO_file );
	 fprintf(fpm,"  deflection exaggeration: %.1f ", exagg );
	 fprintf(fpm,"  load case %d of %d \"\n", lc, nL );

	 fprintf(fpm,"plot '%s' u 2:3 t 'undeformed mesh' w lp ", meshpath);
	 if (!anlyz) fprintf(fpm,"lw 2 lt 1 pt 6 \n");
	 else fprintf(fpm,"lw 1 lt 5 pt 6, '%s' u 1:2 t 'load case %d of %d' w l lw 2 lt 3\n", meshfl, lc, nL );

	 fprintf(fpm,"%c set parametric\n", D3 );
	 fprintf(fpm,"%c set view 60, 70, 1 \n", D3 );
	 fprintf(fpm,"%c set nokey\n", D3 );
	 fprintf(fpm,"%c set xlabel 'x'\n", D3 );
	 fprintf(fpm,"%c set ylabel 'y'\n", D3 );
	 fprintf(fpm,"%c set zlabel 'z'\n", D3 );
/*	 fprintf(fpm,"%c set nolabel\n", D3 );	*/
	 fprintf(fpm,"%c splot '%s' u 2:3:4 t 'undeformed mesh' w lp ",
								D3, meshpath );
	 if (!anlyz) fprintf(fpm," lw 2 lt 1 pt 6 \n");
	 else fprintf(fpm," lw 1 lt 5 pt 6, '%s' u 1:2:3 t 'load case %d of %d' w l lw 2 lt 3\n",meshfl, lc, nL );
	}

	fclose(fpm);

	return;
}


/*------------------------------------------------------------------------------
MODAL_MESH  -  create mesh data of the mode-shape meshes, use gnuplot	19oct98
	 useful gnuplot options: set noxtics noytics noztics noborder view nokey
------------------------------------------------------------------------------*/
void modal_mesh(
		char IO_file[], char meshpath[], char modepath[],
		char plotpath[], char *title,
		int nJ, int nB, int DoF, int nM,
		vec3 *xyz, double *L,
		int *J1, int *J2, float *p,
		double **M, double *f, double **V,
		double exagg, int anlyz
){
	FILE	*fpm;
	double mpfX, mpfY, mpfZ;	/* mode participation factors	*/
	double *msX, *msY, *msZ;
	double *v;		/* a mode-shape vector */

	int	i, j, m,n, X=0, Y=0, Z=0;
	char	D3 = '#', s1[16],  s2[16], modefl[64];


	msX = dvector(1,DoF);
	msY = dvector(1,DoF);
	msZ = dvector(1,DoF);
	v   = dvector(1,DoF);

	for (i=1; i<=DoF; i++) {	/* modal participation factors */
		msX[i] = msY[i] = msZ[i] = 0.0;
		for (j=1; j<=DoF; j+=6) msX[i] += M[i][j];
		for (j=2; j<=DoF; j+=6) msY[i] += M[i][j];
		for (j=3; j<=DoF; j+=6) msZ[i] += M[i][j];
	}

	if (!anlyz) exagg = 0.0;

	for (m=1; m<=nM; m++) {

		strcpy(modefl,modepath);
		s1[0]='-'; s1[1]='\0'; my_itoa(m,s2,2);  strcat(s1,s2);  strcat(modefl,s1);

		if ((fpm = fopen (modefl, "w")) == NULL) {
			printf (" error: cannot open modal mesh file: %s\n", modefl);
			exit(1);
		}

		fprintf(fpm,"# FRAME3DD ANALYSIS RESULTS  http://frame3dd.sf.net/\n");
		fprintf(fpm,"# %s\n", title );
		fprintf(fpm,"# M O D E   S H A P E   D A T A   F O R   M O D E");
		fprintf(fpm,"   %d\t(global coordinates)\n", m );
		fprintf(fpm,"# deflection exaggeration: %.1f\n\n", exagg );
		mpfX = 0.0;	for (i=1; i<=DoF; i++)    mpfX += V[i][m]*msX[i];
		mpfY = 0.0;	for (i=1; i<=DoF; i++)    mpfY += V[i][m]*msY[i];
		mpfZ = 0.0;	for (i=1; i<=DoF; i++)    mpfZ += V[i][m]*msZ[i];
		fprintf(fpm,"# MODE %5d:   f= %lf Hz, T= %lf sec\n", m,f[m],1./f[m]);
		fprintf(fpm,"#\t\tX- modal participation factor = %12.4e \n", mpfX);
		fprintf(fpm,"#\t\tY- modal participation factor = %12.4e \n", mpfY);
		fprintf(fpm,"#\t\tZ- modal participation factor = %12.4e \n", mpfZ);

		for(i=1; i<=DoF; i++)	v[i] = V[i][m];

		fprintf(fpm,"#      X-dsp       Y-dsp       Z-dsp\n\n");

		for(n=1; n<=nB; n++)
			bent_beam ( fpm, J1[n], J2[n], xyz, L[n], p[n], v, exagg );

		for ( j=1; j<=nJ; j++ ) {
			if (xyz[j].x != 0.0) X=1;	/* check for three-dimensional frame */
			if (xyz[j].y != 0.0) Y=1;
			if (xyz[j].z != 0.0) Z=1;
		}

		if ( X && Y && Z ) D3 = ' ';

		fclose(fpm);

		if ((fpm = fopen (plotpath, "a")) == NULL) {
			printf (" error: cannot append plot file: %s\n",plotpath);
			exit(1);
		}
		fprintf(fpm,"pause -1\n");
		fprintf(fpm,"set nolabel\n");
		fprintf(fpm,"set title '%s     mode %d     %lf Hz'\n",IO_file,m,f[m]);
		fprintf(fpm,"plot '%s' u 2:3 t 'undeformed mesh' w l ", meshpath );
		if (!anlyz) fprintf(fpm," lw 2 lt 1 \n");
		else fprintf(fpm," lw 1 lt 5 , '%s' u 1:2 t 'mode-shape %d' w l lw 2 lt 3\n",
								modefl, m );
		fprintf(fpm,"%c pause -1\n", D3 );
		fprintf(fpm,"%c set nokey\n", D3 );
		fprintf(fpm,"%c splot '%s' u 2:3:4 t 'undeformed mesh' w l ",
								D3, meshpath);
		if (!anlyz) fprintf(fpm," lw 2 lt 1 \n");
		else fprintf(fpm," lw 1 lt 5 , '%s' u 1:2:3 t 'mode-shape %d' w l lw 2 lt 3\n",
								modefl, m );

		fclose(fpm);

	}

	free_dvector(msX,1,DoF);
	free_dvector(msY,1,DoF);
	free_dvector(msZ,1,DoF);
	free_dvector(v,1,DoF);
}

/*------------------------------------------------------------------------------
ANIMATE -  create mesh data of animated mode-shape meshes, use gnuplot	16dec98
	 useful gnuplot options: set noxtics noytics noztics noborder view nokey
	 mpeg movie example:   % convert mesh_file-03-f-*.ps mode-03.mpeg
	 ... requires ImageMagick and mpeg2vidcodec packages
------------------------------------------------------------------------------*/
void animate(
	char IO_file[], char meshpath[], char modepath[], char plotpath[],
	char *title,
	int anim[],
	int nJ, int nB, int DoF, int nM,
	vec3 *xyz, double *L, float *p,
	int *J1, int *J2, double *f, double **V,
	double exagg,
	int pan
){
	FILE	*fpm;

	double x_min = 0, x_max = 0,
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
		*v;

	int	fr, i,j, m,n, X=0, Y=0, Z=0, c, CYCLES=3,
		frame_number = 0,
		total_frames;	/* total number of frames in animation */

	char	D3 = '#',
		Movie = '#',	/* use '#' for no-movie  -OR-  ' ' for movie */
		s1[16], s2[16], modefl[64], framefl[64];

	for (j=1; j<=nJ; j++) {		/* check for three-dimensional frame */
		if (xyz[j].x != 0.0) X=1;
		if (xyz[j].y != 0.0) Y=1;
		if (xyz[j].z != 0.0) Z=1;
		if (xyz[j].x < x_min ) x_min = xyz[j].x;
		if (xyz[j].y < y_min ) y_min = xyz[j].y;
		if (xyz[j].z < z_min ) z_min = xyz[j].z;
		if ( x_max < xyz[j].x ) x_max = xyz[j].x;
		if ( y_max < xyz[j].y ) y_max = xyz[j].y;
		if ( z_max < xyz[j].z ) z_max = xyz[j].z;
	}
	if ( X && Y && Z ) D3 = ' ';


	if ((fpm = fopen (plotpath, "a")) == NULL) {
		printf (" error: cannot append plot file: %s\n",plotpath);
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
			x_min-exagg, x_max+exagg );
	   if (y_min != y_max)
		fprintf(fpm,"set yrange [ %lf : %lf ] \n",
	 		y_min-0.2*(y_max-y_min), y_max+0.2*(y_max-y_min) );
	   else fprintf(fpm,"set yrange [ %lf : %lf ] \n",
			y_min-exagg, y_max+exagg );
	   if (z_min != z_max)
	   	fprintf(fpm,"set zrange [ %lf : %lf ] \n",
			z_min-0.2*(z_max-z_min), z_max+0.2*(z_max-z_min) );
	   else fprintf(fpm,"set zrange [ %lf : %lf ] \n",
			z_min-exagg, z_max+exagg );

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

	    strcpy(modefl,modepath);
	    strcpy(framefl,modepath);
	    s1[0] = '-';  s1[1] = '\0';  my_itoa(m,s2,2);  strcat(s1,s2);
	    strcat(framefl,s1);
	    strcat(s1,"."); my_itoa(fr,s2,3); strcat(s1,s2); strcat(modefl,s1);
	    s1[0] = '-'; s1[1] = 'f'; s1[2] = '-'; s1[3] = '\0';
	    my_itoa(frame_number++,s2,3); strcat(s1,s2); strcat(framefl,s1);
	    s1[0] = '.'; s1[1] = 'p'; s1[2] = 's'; s1[3] = '\0';
	    strcat(framefl,s1);

	    if ( D3 == '#' ) {
		fprintf(fpm,"plot '%s' u 2:3 w l lw 1 lt 5, ", meshpath );
	 	fprintf(fpm," '%s' u 1:2 w l lw 2 lt 3 ;", modefl );
	    } else {
	      if (pan)
 	        fprintf(fpm,"%c set view %5.1f, %5.1f, %4.2f \n", D3,
		rot_x_init + (rot_x_final-rot_x_init)*frame_number/total_frames,
		rot_z_init + (rot_z_final-rot_z_init)*frame_number/total_frames,
		zoom_init + (zoom_final-zoom_init)*frame_number/total_frames );
	      fprintf(fpm,"%c splot '%s' u 2:3:4 w l lw 1 lt 5, ",D3,meshpath);
	      fprintf(fpm," '%s' u 1:2:3 w l lw 2 lt 3;", modefl );
	    }
	    if ( fr==0 && c==1 )	fprintf(fpm,"  pause 1.5 \n");
	    else			fprintf(fpm,"  pause 0.05 \n");
	    fprintf(fpm,"%c  load 'saveplot';\n",Movie);
	    fprintf(fpm,"%c  !mv my-plot.ps %s\n", Movie, framefl );
	  }
	  for ( fr = frames-1; fr > 0; fr-- ) {

	    strcpy(modefl,modepath);
	    strcpy(framefl,modepath);
	    s1[0] = '-';  s1[1] = '\0';  my_itoa(m,s2,2);  strcat(s1,s2);
	    strcat(framefl,s1);
	    strcat(s1,"."); my_itoa(fr,s2,3); strcat(s1,s2); strcat(modefl,s1);
	    s1[0] = '-'; s1[1] = 'f'; s1[2] = '-'; s1[3] = '\0';
	    my_itoa(frame_number++,s2,3); strcat(s1,s2); strcat(framefl,s1);
	    s1[0] = '.'; s1[1] = 'p'; s1[2] = 's'; s1[3] = '\0';
	    strcat(framefl,s1);

	    if ( D3 == '#' ) {
	 	fprintf(fpm,"plot '%s' u 2:3 w l lw 1 lt 5, ", meshpath );
		fprintf(fpm," '%s' u 1:2 w l lw 2 lt 3;", modefl );
	    } else {
	      if (pan)
	        fprintf(fpm,"%c set view %5.1f, %5.1f, %4.2f \n", D3,
		rot_x_init + (rot_x_final-rot_x_init)*frame_number/total_frames,
		rot_z_init + (rot_z_final-rot_z_init)*frame_number/total_frames,
		zoom_init + (zoom_final-zoom_init)*frame_number/total_frames );
	      fprintf(fpm,"%c splot '%s' u 2:3:4 w l lw 1 lt 5, ",D3,meshpath);
	      fprintf(fpm," '%s' u 1:2:3 w l lw 2 lt 3;", modefl );
	    }
	    fprintf(fpm,"  pause 0.05 \n");
	    fprintf(fpm,"%c  load 'saveplot';\n",Movie);
	    fprintf(fpm,"%c  !mv my-plot.ps %s\n", Movie, framefl );
	  }
	 }
	 fr = 0;

	 strcpy(modefl,modepath);
	 s1[0] = '-';  s1[1] = '\0';  my_itoa(m,s2,2);  strcat(s1,s2);
	 strcat(s1,".");  my_itoa(fr,s2,3);  strcat(s1,s2);  strcat(modefl,s1);

	 if ( D3 == '#' ) {
	 	fprintf(fpm,"plot '%s' u 2:3 w l lw 2 lt 5, ", meshpath );
		fprintf(fpm," '%s' u 1:2 w l lw 3 lt 3 \n", modefl );
	 } else {
		fprintf(fpm,"%c splot '%s' u 2:3:4 w l lw 2 lt 5, ",D3,meshpath);
		fprintf(fpm," '%s' u 1:2:3 w l lw 3 lt 3 \n", modefl );
	 }

	 i++;
	}
	fclose(fpm);

	v = dvector(1,DoF);

	i = 0;
	while ( (m = anim[i]) != 0 ) {
	  for ( fr=0; fr<=frames; fr++ ) {

	    strcpy(modefl,modepath);
	    s1[0] = '-';  s1[1] = '\0';  my_itoa(m,s2,2);  strcat(s1,s2);
	    strcat(s1,"."); my_itoa(fr,s2,3); strcat(s1,s2); strcat(modefl,s1);

	    if ((fpm = fopen (modefl, "w")) == NULL) {
		printf (" error: cannot open modal mesh file: %s\n", modefl);
		exit(1);
	    }

	    fprintf(fpm,"# FRAME3DD ANALYSIS RESULTS  http://frame3dd.sf.net/\n");
		fprintf(fpm,"# %s\n", title );
	    fprintf(fpm,"# A N I M A T E D   M O D E   S H A P E   D A T A \n");
	    fprintf(fpm,"# deflection exaggeration: %.1f\n", ex );
	    fprintf(fpm,"# MODE %5d: f= %lf Hz  T= %lf sec\n\n",m,f[m],1./f[m]);

	    ex = exagg*cos( PI*fr/frames );

	    for (j=1; j<=DoF; j++)	v[j] = V[j][m];

	    fprintf(fpm,"#      X-dsp       Y-dsp       Z-dsp\n\n");

	    for (n=1; n<=nB; n++)

		bent_beam ( fpm, J1[n], J2[n], xyz, L[n], p[n], v, ex );

	    fclose(fpm);
          }
	  i++;
	}
	free_dvector(v,1,DoF);
	return;
}


/*------------------------------------------------------------------------------
BENT_BEAM  -  computes cubic deflection functions from beam end deflections
and beam end rotations.  Saves deflected shapes to a file.  These bent shapes
are exact for mode-shapes, and for frames loaded at their joints.	22feb99
------------------------------------------------------------------------------*/
void bent_beam(
	FILE *fp, int j1, int j2, vec3 *xyz,
	double L, float p, double *D, double exagg
){
	double	t1, t2, t3, t4, t5, t6, t7, t8, t9, 	/* coord xfmn	*/
		u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12,
		*a, *b, **A,
		s, v, w, dx, dy, dz;
	int	i1, i2, pd;

	A = dmatrix(1,4,1,4);
	a = dvector(1,4);
	b = dvector(1,4);

	coord_trans ( xyz, L, j1, j2,
				&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p );

	i1 = 6*(j1-1);	i2 = 6*(j2-1);

		/* compute beam end deflections in local coordinates */

	u1  = exagg*(t1*D[i1+1] + t2*D[i1+2] + t3*D[i1+3]);
	u2  = exagg*(t4*D[i1+1] + t5*D[i1+2] + t6*D[i1+3]);
	u3  = exagg*(t7*D[i1+1] + t8*D[i1+2] + t9*D[i1+3]);

	u4  = t1*D[i1+4] + t2*D[i1+5] + t3*D[i1+6];
	u5  = t4*D[i1+4] + t5*D[i1+5] + t6*D[i1+6];
	u6  = t7*D[i1+4] + t8*D[i1+5] + t9*D[i1+6];

	u7  = exagg*(t1*D[i2+1] + t2*D[i2+2] + t3*D[i2+3]);
	u8  = exagg*(t4*D[i2+1] + t5*D[i2+2] + t6*D[i2+3]);
	u9  = exagg*(t7*D[i2+1] + t8*D[i2+2] + t9*D[i2+3]);

	u10 = t1*D[i2+4] + t2*D[i2+5] + t3*D[i2+6];
	u11 = t4*D[i2+4] + t5*D[i2+5] + t6*D[i2+6];
	u12 = t7*D[i2+4] + t8*D[i2+5] + t9*D[i2+6];

		/* curve-fitting problem for a cubic polynomial */

	a[1] = u2;		b[1] = u3;
	a[2] = u8;   		b[2] = u9;
	a[3] = exagg*tan(u6);	b[3] = exagg*tan(-u5);
	a[4] = exagg*tan(u12);	b[4] = exagg*tan(-u11);

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
			xyz[j1].x + dx, xyz[j1].y + dy, xyz[j1].z + dz
		);
	}
	fprintf(fp,"\n\n");

	free_dmatrix(A,1,4,1,4);
	free_dvector(a,1,4);
	free_dvector(b,1,4);

	return;
}


/*------------------------------------------------------------------------------
MY_ITOA  -  Convert an integer n to charcters in s, from K&R, 1978,   p. 59-60
... specialized for portability between GNU GCC and DJGPP GCC
------------------------------------------------------------------------------*/
void my_itoa(int n, char s[], int k){
	int	c, i, j, sign;

	if ((sign = n) < 0) 		/* record sign */
		n = -n;			/* make n positive */
	i = 0;
	do {				/* generate digits in reverse order */
		s[i++] = n % 10 + '0';	/* get next digit */
	} while ((n /= 10) > 0);	/* delete it */
	for (;i<k;)	s[i++] = '0';	/* add leading '0' */
	if (sign < 0)
		s[i++] = '-';
	s[i] = '\0';
					/* reverse order of string s */
	j = 0;
	while ( s[j] != '\0' )	j++;	/* j is length of s - 1 */
	--j;

	for (i = 0; i < j; i++, j--) {
		c = s[i];
		s[i] = s[j];
		s[j] = c;
	}
	return;
}


/*------------------------------------------------------------------------------
GET_FILE_EXT  -  get the file extension,
		return 1 if the extension is ".csv"
		return 2 if the extension is ".fmm"
		return 0 otherwise
------------------------------------------------------------------------------*/
int get_file_ext( char *filename, char *ext )
{
	int	i=0, full_len=0, len=0;

	while ( filename[len++] != '\0' ) /* the length of file filename */ ;
        full_len = len;
	while ( filename[len--] != '.' ) /* the last '.' in filename */ ;
	++len;

	for ( i=0; len < full_len; i++,len++ ) ext[i] = tolower(filename[len]);

	/* test */
//	printf(" filename '%s' has length %d and extension = '%s' \n",
//							filename, len, ext);
//	printf(" Is .CSV? ... = %d \n", !strcmp(ext,".csv") );

	if ( !strcmp(ext,".csv") ) return (1);
	if ( !strcmp(ext,".fmm") ) return (2);
	return(0);
}


/*------------------------------------------------------------------------------
DOTS  -  print a set of dots (periods)
------------------------------------------------------------------------------*/
void dots( int n ){
	int i;
	for (i=1; i<=n; i++)	printf(".");
}

