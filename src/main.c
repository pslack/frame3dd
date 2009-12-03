/*
 FRAME3DD:
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
*//** @file
	Main FRAME3DD program driver
*//** @mainpage
FRAME3DD: a program for static and dynamic structural analysis of 2D and 3D
frames and trusses with elastic and geometric stiffness.

Also included is a system for parsing Microstran .arc 'Archive' files and
for parsing calculated force and displacement output files (.p1 format) from
Microstran. See @ref mstranp. It is intended that ultimately the .arc format
be an alternative method of inputting data to the FRAME3DD program,
but currently these two parts of the code are distinct.

For more information go to http://frame3dd.sourceforge.net/

The input file format for FRAME is defined in doc/user_manual.html

Henri P. Gavin hpgavin@duke.edu (main FRAME3DD code) <br>
John Pye john.pye@anu.edu.au (Microstran parser and viewer)

For compilation/installation, see README.txt.

*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "frame3dd.h"
#include "frame3dd_io.h"
#include "eig.h"
#include "nrutil.h"


// compile the Frame3DD analysis into another code, such as a GUI
#ifdef WITH_GLOBALS	
 int run_kernel ( int argc, char *argv[] ) {
#endif

// compile Frame3DD to run as a stand-alone code through the terminal
#ifndef WITH_GLOBALS
 int main ( int argc, char *argv[] ) {
#endif

	char	IN_file[FILENMAX],	/* the input  data filename	*/
		OUT_file[FILENMAX],	/* the output data filename	*/
		title[256],		/* the title of the analysis	*/
		meshpath[FRAME3DD_PATHMAX] = "EMPTY_MESH", /* mesh data path */
		plotpath[FRAME3DD_PATHMAX] = "EMPTY_PLOT", /* plot file path */
		modepath[FRAME3DD_PATHMAX] = "EMPTY_MODE", /* mode data path */
		temppath[FRAME3DD_PATHMAX] = "EMPTY_TEMP"; /* temp data path */

	FILE	*fp;		/* input and output file pointer	*/

	vec3	*xyz;		/* X,Y,Z joint coordinates (global)	*/

	float	*r,		/* joint size radius, for finite sizes	*/
		*Ax,*Asy, *Asz,	/* cross section areas, incl. shear	*/
		*J,*Iy,*Iz,	/* section inertias			*/
		*E, *G,		/* elastic modulus and shear moduli	*/
		*p,		/* roll of each member, radians		*/
		***U,		/* uniform distributed member loads	*/
		***W,		/* trapizoidal distributed member loads	*/
		***P,		/* member concentrated loads		*/
		***T,		/* member temperature  loads		*/
		**Dp,		/* prescribed joint displacements	*/
		*d, *BMs,	/* member densities and extra inertia	*/
		*JMs, 		/* mass of a joint			*/
		*JMx,*JMy,*JMz,	/* inertia of a joint in global coord	*/
		pan=1.0;	/* >0: pan during animation; 0: don't	*/

	double	**K, **Ks=NULL,	/* global stiffness matrix		*/
		traceK = 0.0,	/* trace of the global stiffness matrix	*/
		**M = NULL,	/* global mass matrix			*/
		traceM = 0.0,	/* trace of the global mass matrix	*/
		**Fo_mech,	/* mechanical load vectors,  load cases	*/
		**Fo_temp,	/* thermal load vectors, all load cases	*/
		**Fo, *F,	/* general load vectors			*/
		***feF_mech,	/* fixed end forces from mech loads	*/
		***feF_temp,	/* fixed end forces from temp loads	*/
		**feF,		/* a general set of fixed end forces	*/
		*D, *dD,	/* displacement and displ increment	*/
		//dDdD = 0.0,	/* dD' * dD				*/
		*Fe = NULL,	/* equilibrium error in nonlinear anlys	*/
		*L  = NULL,	/* joint-to-joint length of each element*/
		*Le = NULL,	/* effcve lngth, accounts for joint size*/
		**Q = NULL,	/* local member joint end-forces	*/
		tol = 1.0e-9,	/* tolerance for modal convergence	*/
		shift = 0.0,	/* shift-factor for rigid-body-modes	*/
		struct_mass,	/* mass of structural system		*/
		total_mass,	/* total structural mass and extra mass */
		*f  = NULL,	/* resonant frequencies			*/
		**V = NULL,	/* resonant mode-shapes			*/
		error = 1.0,	/* rms equilibrium error and reactions	*/
		Cfreq = 0.0,	/* frequency used for Guyan condensation*/
		**Kc, **Mc,	/* condensed stiffness and mass matrices*/
		exagg_static=10,/* exaggerate static displ. in mesh data*/
		exagg_modal=10;	/* exaggerate modal displ. in mesh data	*/

	int	nJ=0,		/* number of Joints 			*/
		nE=0,		/* number of frame Elements		*/
		nL=0, lc=0,	/* number of Load cases			*/
		DoF=0, i, j, n,	/* number of Degrees of Freedom		*/
		nR=0,		/* number of restrained joints		*/
		nD[_NL_],	/* number of prescribed nodal displ'nts	*/
		nF[_NL_],	/* number of loaded joints		*/
		nU[_NL_],	/* number of members w/ unifm dist loads*/
		nW[_NL_],	/* number of members w/ trapz dist loads*/
		nP[_NL_],	/* number of members w/ conc point loads*/
		nT[_NL_],	/* number of members w/ temp. changes	*/
		nI=0,		/* number of joints w/ extra inertia	*/
		nC=0,		/* number of condensed joints		*/
		*J1, *J2,	/* begin and end joint numbers		*/
		shear=0,	/* indicates shear deformationi		*/
		geom=0,		/* indicates  geometric nonlinearity	*/
		anlyz=1,	/* 1: stiffness analysis, 0: data check	*/
		*R, sumR,	/* reaction data, total no. of reactions*/
		nM=0,		/* number of desired modes		*/
		Mmethod,	/* 1: Subspace Jacobi, 2: Stodola	*/
		nM_calc,	/* number of modes to calculate		*/
		lump=1,		/* 1: lumped, 0: consistent mass matrix */
		iter=0,		/* number of iterations			*/
		ok=1,		/* number of (-ve) diag. terms of L D L' */
		anim[20],	/* the modes to be animated		*/
		Cdof=0,		/* number of condensed degrees o freedom*/
		Cmethod=0,	/* matrix condensation method		*/
		*q,		/* vector of DoF's to condense		*/
		*m,		/* vector of modes to condense		*/
		filetype=0,	/* 1 if .CSV, 2 if file is Matlab	*/
		debug=0,	/* 1: debugging screen output, 0: none	*/
		verbose=1;	/* 1: copious screen output, 0: none	*/

	int	shear_flag= -1,	/*   over-ride input file value		*/
		geom_flag = -1,	/*   over-ride input file value		*/
		anlyz_flag= -1,	/*   over-ride input file value		*/
		lump_flag = -1,	/*   over-ride input file value		*/
		modal_flag= -1,	/*   over-ride input file value		*/
		write_matrix=-1,/*   write stiffness and mass matrix	*/
		condense_flag=-1; /* over-ride input file value		*/

	int	sfrv=0;		/* *scanf return value for err checking	*/

	double	exagg_flag=-1.0, /*  over-ride input file value		*/
		tol_flag  =-1.0, /*  over-ride input file value		*/
		shift_flag=-1.0; /*  over-ride input file value		*/

	float	pan_flag = -1.0; /*  over-ride input file value		*/

	char	extn[16];	/* Input Data file name extension	*/


	parse_options ( argc, argv, IN_file, OUT_file, 
			&shear_flag, &geom_flag, &anlyz_flag, &exagg_flag, 
			&lump_flag, &modal_flag, &tol_flag, &shift_flag, 
			&pan_flag, &write_matrix, &condense_flag, &verbose, &debug);

	if ( verbose ) {
		fprintf(stderr,"\n FRAME3DD version: %s\n", VERSION);
		fprintf(stderr," Analysis of 2D and 3D structural frames with elastic and geometric stiffness.\n");
		fprintf(stderr," http://frame3dd.sf.net\n");
		fprintf(stderr," GPL Copyright (C) 1992-2009, Henri P. Gavin\n");
		fprintf(stderr," This is free software with absolutely no warranty.\n");
		fprintf(stderr," For details, see the GPL license file, LICENSE.txt\n\n");
	}

	/* open the input data file */

	if ((fp = fopen (IN_file, "r")) == NULL) {
		fprintf (stderr,"\n ERROR: cannot open file '%s'\n\n", IN_file);
		display_help();
		if ( argc == 1 ) {
			fprintf(stderr," Press the 'Enter' key to close.\n");
			(void) getchar();	// clear the buffer ?? 
			while( !getchar() ) ;	// wait for the Enter key 
		}
		exit(1);
	}

	filetype = get_file_ext( IN_file, extn ); /* .CSV or .FMM or other? */

//	temp_file_location("frame3dd.3dd",temppath,FRAME3DD_PATHMAX);
	output_path("frame3dd.3dd",temppath,FRAME3DD_PATHMAX,NULL);

	parse_input(fp, temppath);	/* strip comments from input data */
	fclose(fp);

	/* open the clean input file */
	if ((fp = fopen (temppath, "r")) == NULL) {
		fprintf (stderr,"\n ERROR: cannot open cleaned input file '%s'\n",temppath);
		exit(1);
	}

	frame3dd_getline(fp, title, 256);
	if (verbose ) fprintf(stderr," ** %s ** \n\n", title );

	sfrv=fscanf(fp, "%d", &nJ );		/* number of joints	*/
	if (sfrv != 1)	sferr("nJ value for number of joints");
	if ( verbose ) {
	 printf(" number of joints "); dots(stdout,35); printf(" nJ = %3d ",nJ);
	}

					/* allocate memory for joint data ... */
	r   =  vector(1,nJ);		/* rigid radius around each joint */
	xyz = (vec3 *)malloc(sizeof(vec3)*(1+nJ));	/* joint coordinates */

	read_joint_data ( fp, nJ, xyz, r );
	if ( verbose )	printf(" ... complete\n");

	DoF = 6*nJ;		/* total number of degrees of freedom	*/

	R   = ivector(1,DoF);	/* allocate memory for reaction data ... */
	read_reaction_data ( fp, DoF, nJ, &nR, R, &sumR, verbose );
	if ( verbose )	printf(" ... complete\n");

	sfrv=fscanf(fp, "%d", &nE );	/* number of frame elements	*/
	if (sfrv != 1)	sferr("nE value for number of frame elements");
	if ( verbose ) {
	 printf(" number of frame elements"); dots(stdout,28); printf(" nE = %3d ",nE);
	}
	if ( nJ > nE + 1) {
		fprintf(stderr,"warning: %d joints and %d members...", nJ, nE );
		fprintf(stderr," not enough members to connect all joints.\n");
    	}

				/* allocate memory for frame elements ... */
	L   = dvector(1,nE);	/* length of each element		*/
	Le  = dvector(1,nE);	/* effective length of each element	*/

	J1  = ivector(1,nE);	/* joint #1 of each element		*/
	J2  = ivector(1,nE);	/* joint #2 of each element		*/

	Ax  =  vector(1,nE);	/* cross section area of each element	*/
	Asy =  vector(1,nE);	/* shear area in local y direction 	*/
	Asz =  vector(1,nE);	/* shear area in local z direction	*/
	J   =  vector(1,nE);	/* torsional moment of inertia 		*/
	Iy  =  vector(1,nE);	/* bending moment of inertia about y-axis */
	Iz  =  vector(1,nE);	/* bending moment of inertia about z-axis */

	E   =  vector(1,nE);	/* frame element Young's modulus	*/
	G   =  vector(1,nE);	/* frame element shear modulus		*/
	p   =  vector(1,nE);	/* member rotation angle about local x axis */

	read_frame_element_data( fp, nJ, nE, xyz,r, L, Le, J1, J2,
					Ax, Asy, Asz, J, Iy, Iz, E, G, p );
	if ( verbose) 	printf(" ... complete\n");


	read_run_data ( fp, OUT_file, &shear, shear_flag, &geom, geom_flag,
			meshpath, plotpath,
			&exagg_static, exagg_flag, &anlyz, anlyz_flag, debug );

	sfrv=fscanf(fp, "%d", &nL );	/* number of load cases		*/
	if (sfrv != 1)	sferr("nL value for number of load cases");
	if ( verbose ) {
		printf(" number of load cases ");
		dots(stdout,31); printf(" nL = %3d \n",nL);
	}

	if ( nL < 1 ) {
		fprintf(stderr,"\n ERROR: the number of load cases must be at least 1\n");
		exit(1);
	}
	if ( nL >= _NL_ ) {
		fprintf(stderr,"\n ERROR: maximum of %d load cases allowed\n", _NL_-1);
		exit(1);
	}
					/* allocate memory for loads ... */
	U   =  D3matrix(1,nL,1,nE,1,4);    /* uniform load on each member */
	W   =  D3matrix(1,nL,1,10*nE,1,13);/* trapezoidal load on each member */
	P   =  D3matrix(1,nL,1,10*nE,1,5); /* internal point load each member */
	T   =  D3matrix(1,nL,1,nE,1,8);    /* internal temp change each member*/
	Dp  =  matrix(1,nL,1,DoF); /* prescribed displacement of each joint */

	Fo_mech  = dmatrix(1,nL,1,DoF);	/* mechanical load vector	*/
	Fo_temp  = dmatrix(1,nL,1,DoF);	/* temperature load vector	*/
	Fo  = dmatrix(1,nL,1,DoF);	/* external load vector		*/
	F   = dvector(1,DoF);	/* external load vector with react'ns	*/

	feF_mech =  D3dmatrix(1,nL,1,nE,1,12); /* feF due to mech loads */
	feF_temp =  D3dmatrix(1,nL,1,nE,1,12); /* feF due to temp loads */
	feF      = dmatrix(1,nE,1,12);	/* fixed end forces		*/

	K   = dmatrix(1,DoF,1,DoF);	/* global stiffness matrix	*/
	Q   = dmatrix(1,nE,1,12);	/* end forces for each member	*/

	D   = dvector(1,DoF);	/* displacments of each joint		*/
	dD  = dvector(1,DoF);	/* incremental displ. of each joint	*/

	d   =  vector(1,nE);	/* mass density for each member		*/
	BMs =  vector(1,nE);	/* lumped mass for each frame element	*/
	JMs =  vector(1,nJ);	/* joint mass for each joint		*/
	JMx =  vector(1,nJ);	/* joint inertia about global X axis	*/
	JMy =  vector(1,nJ);	/* joint inertia about global Y axis	*/
	JMz =  vector(1,nJ);	/* joint inertia about global Z axis	*/

	q = ivector(1,DoF); 	/* vector of condensed degrees of freedom */
	m = ivector(1,DoF); 	/* vector of condensed mode numbers	*/


	read_and_assemble_loads( fp, nJ, nE, nL, DoF, xyz, L, Le, J1, J2,
				Ax,Asy,Asz, Iy,Iz, E, G, p, R, shear,
				nF, nU, nW, nP, nT, nD,
				Q, Fo_mech, Fo_temp, U, W, P, T,
				Dp, feF_mech, feF_temp, verbose );

	for (i=1; i<=DoF; i++){
		for (lc=1; lc<=nL; lc++){
			Fo[lc][i] = Fo_temp[lc][i] + Fo_mech[lc][i];
		}
	}
	if ( verbose ) {
		printf("                                                     ");
		printf(" load data ... complete\n");
	}

	read_mass_data( fp, IN_file, nJ, nE, &nI, d, BMs, JMs, JMx, JMy, JMz,
			L, Ax, &total_mass, &struct_mass, &nM,
			&Mmethod, modal_flag, 
			&lump, lump_flag, &tol, tol_flag, &shift, shift_flag,
			&exagg_modal,
			modepath,
			anim, &pan, pan_flag, 
			verbose, debug );
	if ( verbose ) {
		printf("                                                     ");
		printf(" mass data ... complete\n");
	}

	read_condensation_data( fp, nJ,nM, &nC, &Cdof, 
			&Cmethod, condense_flag, q,m, verbose );

	if( nC>0 && verbose ) {
		printf("                                      ");
		printf(" matrix condensation data ... complete\n");
	}

	fclose(fp);

	/* open the output data file for appending */

	fp = fopen(OUT_file, "a");
	if(fp==NULL){
		fprintf(stderr,"Unable to append to output file '%s'!\n",
								OUT_file);
		exit(1);
	}

	write_input_data ( fp, title, nJ,nE,nL, nD,nR, nF,nU,nW,nP,nT,
				xyz, r, J1,J2, Ax,Asy,Asz, J,Iy,Iz, E,G, p,
				Fo, Dp, R, U, W, P, T, shear, anlyz, geom );

	if ( anlyz ) {			/* solve the problem	*/
		srand(time(NULL));
		for (lc=1; lc<=nL; lc++) {	/* begin load case loop	*/

			if ( verbose )
				printf("\n Load Case %d of %d ... \n", lc,nL );

			for (i=1; i<=DoF; i++)	D[i] = dD[i] = 0.0;

			assemble_K ( K, DoF, nE, xyz,r, L, Le, J1, J2,
				Ax, Asy, Asz, J,Iy,Iz, E, G, p, shear, geom, Q );

#ifdef MATRIX_DEBUG
			save_dmatrix ( DoF, DoF, K, "Kf" ); /* free stiffness matrix */
#endif

			/* apply temperature loads first ... */
			if (nT[lc] > 0) {
				if ( verbose )
					printf(" Linear Elastic Analysis ... Temperature Loads\n");
				apply_reactions ( DoF, R, Dp[lc], Fo_temp[lc], F, K, 't' );
				solve_system( K, dD, F, DoF, &ok, verbose );
				for (i=1; i<=DoF; i++)	D[i] += dD[i];
				end_forces ( Q, nE, xyz, L, Le, J1,J2,
					Ax, Asy,Asz, J,Iy,Iz, E,G, p, D, shear, geom );
			}

			assemble_K ( K, DoF, nE, xyz, r, L, Le, J1, J2,
				Ax, Asy, Asz, J,Iy,Iz, E, G, p, shear, geom, Q );

			/* then add mechanical loads ... */
			if ( nF[lc] > 0 || nU[lc] > 0 || nW[lc] > 0 || nP[lc] > 0 || nD[lc] > 0 ) {
				if ( verbose )
					printf(" Linear Elastic Analysis ... Mechanical Loads\n");
				apply_reactions ( DoF, R, Dp[lc], Fo_mech[lc], F, K, 'm' );
				solve_system( K, dD, F, DoF, &ok, verbose );
				for (i=1; i<=DoF; i++)	D[i] += dD[i];
				end_forces ( Q, nE, xyz, L, Le, J1,J2,
					Ax, Asy,Asz, J,Iy,Iz, E,G, p, D, shear, geom );
			}

			/* Newton-Raphson iterations for geometric non-linearity */
			if ( geom ) {
				Fe  = dvector( 1, DoF );	/* force equilibrium error  */
				Ks  = dmatrix( 1, DoF, 1, DoF );
				for (i=1;i<=DoF;i++){ /* initialize Broyden secant stiffness */
					for(j=i;j<=DoF;j++){
						Ks[i][j]=Ks[j][i]=K[i][j];
					}
				}
				if ( verbose )
					printf("\n Non-Linear Elastic Analysis ...\n");
			}

			/* Newton Raphson iteration for geometric nonlinearity */
			ok = 0; iter = 0; error = 1.0;	/* re-initialize */
			while ( geom && error > tol && iter < 10 && ok >= 0) {
				++iter;

				assemble_K ( K, DoF, nE, xyz, r, L, Le, J1, J2,
					Ax, Asy, Asz, J,Iy,Iz, E, G, p, shear, geom, Q );

				apply_reactions ( DoF, R, Dp[lc], Fo[lc], F, K, 'm' );

				for (i=1; i<=DoF; i++) {	/* equilibrium error	*/
					Fe[i] = F[i];
					for (j=1; j<=DoF; j++)
						if ( K[i][j] != 0.0 && D[j] != 0.0 )
							Fe[i] -= K[i][j]*D[j];
				}

				/* Broyden secant stiffness matrix update */
/*
				dDdD = 0.0;
				for (i=1; i<=DoF; i++) dDdD += dD[i]*dD[i];
				for (i=1; i<=DoF; i++)
					for (j=1; j<=DoF; j++)
						Ks[i][j] -= Fe[i]*dD[j] / dDdD;
*/

				apply_reactions ( DoF, R, Dp[lc], Fe, Fe, Ks, 'm' );

				solve_system ( K, dD, Fe, DoF, &ok, verbose );

				if ( ok < 0 ) {
					fprintf(stderr,"   The stiffness matrix is not pos-def. \n");
					fprintf(stderr,"   Reduce loads and re-run the analysis.\n");
					break;
				}			/* end Newton Raphson iterations */

				for (i=1; i<=DoF; i++)	D[i] += dD[i];	/* increment D */

				end_forces ( Q, nE, xyz, L, Le,
					J1,J2, Ax, Asy,Asz, J,Iy,Iz, E,G, p, D, shear, geom );

						 /* convergence criteria:  */
				//error = rel_norm ( dD, D, DoF ); /* displacement increment */
				error = rel_norm ( Fe, F, DoF ); /* force balance	   */

				if ( verbose ) { 
					printf("   NR iteration %3d ---", iter);
					printf(" RMS equilibrium precision: %8.2e \n", error);
				}
			}

			if ( geom ) {
				free_dvector(Fe, 1, DoF );
				free_dmatrix(Ks, 1, DoF, 1, DoF );
			}

			if ( write_matrix ) /* write static stiffness matrix */
				save_ut_dmatrix ( DoF, K, "Ks" );

			for (i=1; i<=12; i++)
			    for (n=1; n<=nE; n++)
				feF[n][i] = feF_temp[lc][n][i] + feF_mech[lc][n][i];

			equilibrium ( xyz, L, J1,J2, Fo[lc], R, p, Q, feF,
						nE, DoF, &error, verbose );

			if ( verbose ) {
				printf("  RMS relative equilibrium precision: %9.3e", error );
				evaluate ( error );
			}


			write_static_results ( fp, nJ,nE,nL,lc, DoF, J1,J2, Fo[lc],
									 D,R,Q, error, ok );

			if ( filetype == 1 ) {		// .CSV format output
				write_static_csv(OUT_file, title,
					nJ,nE,nL,lc, DoF, J1,J2, Fo[lc], D,R,Q, error, ok );
			}

			if ( filetype == 2 ) {		// matlab format output
				write_static_mfile (OUT_file, title,
					nJ,nE,nL,lc, DoF, J1,J2, Fo[lc], D,R,Q, error, ok );
			}

			if ( verbose ) 
				printf("\n    If the program pauses here for very long,"
				" hit CTRL-C to stop execution, \n"
				"    reduce exagg_static in the Input Data,"
				" and re-run the analysis. \n");

			static_mesh ( IN_file, meshpath, plotpath, title, nJ, nE, nL, lc,
				DoF, xyz, L, J1,J2, p, D, exagg_static, anlyz);

		} /* end load case loop */
	} else {
	
		if ( verbose ) {
			printf("\n * %s *\n", title );
			printf("  DATA CHECK ONLY.\n");
		}
		static_mesh ( IN_file, meshpath, plotpath, title, nJ, nE, nL, lc,
				DoF, xyz, L, J1,J2, p, D, exagg_static, anlyz);
	}


	if(nM > 0){ /* modal analysis */

		if ( verbose ) printf("\n Modal Analysis ...\n");

		nM_calc = (nM+8)<(2*nM) ? nM+8 : 2*nM;

		M   = dmatrix(1,DoF,1,DoF);
		f   = dvector(1,nM_calc);
		V   = dmatrix(1,DoF,1,nM_calc);

		assemble_M ( M, DoF, nJ, nE, xyz, r, L, J1,J2,
				Ax, J,Iy,Iz, p, d, BMs, JMs, JMx, JMy, JMz, lump );

#ifdef MATRIX_DEBUG
		save_dmatrix ( DoF, DoF, M, "Mf" );	/* free mass matrix */
#endif

		for (j=1; j<=DoF; j++) {
			if ( !R[j] ) {
				traceK += K[j][j];
				traceM += M[j][j];
			}
		}
		for (i=1; i<=DoF; i++) {
			if ( R[i] ) {	/* apply reactions to upper triangle */
				K[i][i] = traceK * 1e2;
				M[i][i] = traceM;
				for (j=i+1; j<=DoF; j++)
					K[j][i]=K[i][j]=M[j][i]=M[i][j] = 0.0;
		    }
		}

		if ( write_matrix ) {
			save_ut_dmatrix ( DoF, K, "Kd" );/* dynamic stff matx */
			save_ut_dmatrix ( DoF, M, "Md" );/* dynamic mass matx */
		}

		if(anlyz) {
			if( Mmethod == 1 )
				subspace( K, M, DoF, nM_calc, f, V, tol,shift,&iter,&ok, verbose );
			if( Mmethod == 2 )
				stodola ( K, M, DoF, nM_calc, f, V, tol,shift,&iter,&ok, verbose );

			for (j=1; j<=nM_calc; j++)	f[j] = sqrt(f[j])/(2.*PI);

			write_modal_results ( fp, nJ,nE,nI, DoF, M,f,V,
					total_mass, struct_mass,
					iter, sumR, nM, shift, lump, tol, ok );
		}
	}

	fprintf(fp,"\n");
	fclose (fp);

	if(nM > 0 && anlyz){

		modal_mesh ( IN_file, meshpath, modepath, plotpath, title,
		    nJ,nE, DoF, nM, xyz, L, J1,J2, p, M,f,V,exagg_modal,anlyz);
		animate ( IN_file, meshpath, modepath, plotpath, title,anim,
		    nJ,nE, DoF, nM, xyz, L, p, J1,J2, f,V, exagg_modal, pan );
	}

	if(nC > 0){		/* matrix condensation of stiffness and mass */

		if ( verbose ) printf("\n Matrix Condensation ...\n");

		if(Cdof > nM && Cmethod == 3){
			fprintf(stderr,"  Cdof > nM ... Cdof = %d  nM = %d \n",
				 Cdof, nM );
			fprintf(stderr,"  The number of condensed degrees of freedom");
			fprintf(stderr," may not exceed the number of computed modes");
			fprintf(stderr," when using dynamic condensation.\n");
			exit(1);
		}

		Kc = dmatrix(1,Cdof,1,Cdof);
		Mc = dmatrix(1,Cdof,1,Cdof);

		if ( m[1] > 0 && nM > 0 )	Cfreq = f[m[1]];

		if ( Cmethod == 1 ) {	/* static condensation only	*/
			condense(K, DoF, q, Cdof, Kc, verbose );
			if ( verbose )
				printf("   static condensation of K complete\n");
		}
		if ( Cmethod == 2 ) {
			guyan(M, K, DoF, q, Cdof, Mc,Kc, Cfreq, verbose );
			if ( verbose ) {
				printf("   Guyan condensation of K and M complete");
				printf(" ... dynamics matched at %f Hz.\n", Cfreq );
			}
		}
		if ( Cmethod == 3 && nM > 0 ) {
			dyn_conden(M,K, DoF, R, q, Cdof, Mc,Kc, V,f, m, verbose );
			if ( verbose ) 
				printf("   dynamic condensation of K and M complete\n");
		}
		save_dmatrix(Cdof, Cdof, Kc, "Kc" );
		save_dmatrix(Cdof, Cdof, Mc, "Mc" );

		free_dmatrix(Kc, 1,Cdof,1,Cdof );
		free_dmatrix(Mc, 1,Cdof,1,Cdof );
	}

	/* deallocate memory used for each frame analysis variable */
	deallocate ( nJ, nE, nL, nF, nU, nW, nP, nT, DoF, nM,
			xyz, r, L, Le, J1, J2, R,
			Ax, Asy, Asz, J, Iy, Iz, E, G, p,
			U,W,P,T, Dp, Fo_mech, Fo_temp,
			feF_mech, feF_temp, feF, Fo, F,
			K, Q, D, dD,
			d,BMs,JMs,JMx,JMy,JMz, M,f,V, q, m
	);

	if ( verbose ) printf("\n");

	/* wait for keyboard entry to close the terminal */
	if ( argc == 1 ) {
	 fprintf(stderr," The Output Data was appended to %s \n", OUT_file );
	 fprintf(stderr," A Gnuplot script was written to %s \n", plotpath );
	 fprintf(stderr," Press the 'Enter' key to close.\n");
	 (void) getchar();	// clear the buffer ?? 
	 while( !getchar() ) ;	// wait for the Enter key to be hit 
	}

	return(0);
}
