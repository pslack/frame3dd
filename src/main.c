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
*//** @file
	Main FRAME3DD program driver
*//** @mainpage
	FRAME3DD: a program for static and dynamic structural analysis of 2D and 3D
	frames and trusses with elastic and geometric stiffness.

	Also included is a system for parsing Microstran .arc 'Archive' files and
	for parsing calculated force and displacement output files (.p1 format) from
	Microstran. See @ref mstranp. It is intended that ultimately the .arc format be an alternative
	method of inputting data to the FRAME3DD program, but currently these two
	parts of the code are distinct.

	For more information go to http://frame3dd.sourceforge.net/

	The input file format for FRAME is defined in doc/user_manual.html

	Henri P. Gavin hpgavin@duke.edu (main FRAME3DD code) <br>
	John Pye john.pye@anu.edu.au (Microstran parser and viewer)

	For compilation/installation, see README.txt.
*/
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "nrutil.h"

/* must come after the above, because of the sneaky #defines in common.h */
#include "frame3dd.h"
#include "common.h" 
#include "coordtrans.h"
#include "ldl_dcmp.h"
#include "eig.h"
#include "frame3dd_io.h"

#define _NL_ 128

int main(int argc, char *argv[]){

	char	IO_file[96],	/* the input/output filename		*/
		title[256],	/* the title of the analysis		*/
		mesh_file[96],	/* frame mesh data filename		*/
		plot_file[96],	/* frame mesh plot filename		*/
		mode_file[96];	/* mode-shape mesh data filename	*/

	FILE	*fp;		/* input/output file pointer		*/

	vec3	*pos;		/* X,Y,Z joint coordinates (global)	*/

	double *r,		/* joint size radius, for finite sizes	*/
		**K, **Ks,	/* global stiffness matrix		*/
		traceK = 0.0,	/* trace of the global stiffness matrix	*/
		**M = NULL,	/* global mass matrix			*/
		traceM = 0.0,	/* trace of the global mass matrix	*/
		**Fo_mech,	/* mechanical load vector		*/
		**Fo_temp,	/* thermal load vector		*/
		*Fo_mech_lc,	/* mechanical load vector		*/
		*Fo_temp_lc,	/* thermal load vector		*/
		**Fo, *F_lc, *Fo_lc,	/* general load vectors		*/
		***feF_mech,	/* fixed end forces from mech loads	*/
		***feF_temp,	/* fixed end forces from temp loads	*/
		**feF,		/* a general set of fixed end forces	*/
		*D, *dD,	/* displacement and displ increment	*/
		/*dDdD = 0.0,*/	/* dD' * dD				*/
		*Fe,		/* equilibrium error in nonlinear anlys	*/
		*Dp,		/* prescribed joint displacements	*/
		***W,		/* uniform distributed member loads	*/
		***P,		/* member concentrated loads		*/
		***T,		/* member temperature  loads		*/
		*L, *E, *G,	/* length, elastic and shear modulous	*/
		*p,		/* roll of each member, radians		*/
		*Le,		/* effcve lngth, accounts for joint size*/
		*Ax,*J,*Iy,*Iz,	/* area and inertias member coordinates	*/
		*Asy, *Asz,	/* shear areas for shear deformations	*/
		**Q,		/* local member joint end-forces	*/
		*d, *BMs,	/* member densities and extra inertia	*/
		*JMs, 		/* mass of a joint			*/
		struct_mass,	/* mass of structural system		*/
		total_mass,	/* total structural mass and extra mass */
		*JMx,*JMy,*JMz,	/* inertia of a joint in global coord	*/
		tol = 1e-5,	/* tolerance for modal convergence	*/
		shift		/* shift-factor for rigid-body-modes	*/
		,*f = NULL  /* resonant frequencies */
		,**V = NULL,/* resonant mode-shapes	*/
		error = 1.0,	/* rms equilibrium error and reactions	*/
		Cfreq = 0.0,	/* frequency used for Guyan condensation*/
		**Kc, **Mc,	/* condensed stiffness and mass matrices*/
		exagg;		/* exaggerate deformations in mesh data	*/

	int	nJ,		/* number of Joints 			*/
		nM,		/* number of Members			*/
		nL, lc,		/* number of Load cases			*/
		DoF, i, j, n,	/* number of Degrees of Freedom		*/
		nR,		/* number of restrained joints		*/
		nD,		/* number of prescribed nodal displ'nts	*/
		nF[_NL_],	/* number of loaded joints		*/
		nW[_NL_],	/* number of members w/ unif distr loads*/
		nP[_NL_],	/* number of members w/ conc point loads*/
		nT[_NL_],	/* number of members w/ temp. changes	*/
		nI,             /* number of joints w/ extra inertia	*/
		nC,		/* number of condensed joints		*/
		*J1, *J2,	/* begin and end joint numbers		*/
		shear,		/* indicates shear deformationi		*/
		geom=0,		/* indicates  geometric nonlinearity	*/
		anlyz=1,	/* 1: stiffness analysis, 0: data check	*/
		*R, sumR,	/* reaction data, total no. of reactions*/
		modes,		/* number of desired modes		*/
		Mmethod,	/* 1: Subspace Jacobi, 2: Stodola	*/
		calc_modes,	/* number of modes to calculate		*/
		lump=1,		/* 1: lumped, 0: consistent mass matrix */
		iter=0,		/* number of iterations			*/
		ok=1,		/* number of (-ve) diag. terms of L D L' */
		anim[20],	/* the modes to be animated		*/
		pan=1,		/* 1: pan during animation; 0: don't	*/
		Cdof,		/* number of condensed degrees o freedom*/
		Cmethod,	/* matrix condensation method		*/
		*q,		/* vector of DoF's to condense		*/
		*m;		/* vector of modes to condense		*/

    fprintf(stderr," FRAME3DD version: %s\n", VERSION);
    fprintf(stderr," GPL Copyright (C) 1992-2008, Henri P. Gavin\n");
    fprintf(stderr," http://frame3dd.sf.net\n");
    fprintf(stderr," This is free software with absolutely no warranty.\n");
    fprintf(stderr," For details, see LICENSE.txt\n\n");

	if (argc < 2) {
		fprintf (stderr," Please enter the input/output file name: ");
		scanf("%s", IO_file );
		fprintf (stderr," You entered file name: %s \n", IO_file );
	} else strcpy( IO_file , argv[1] );

	if ((fp = fopen (IO_file, "r")) == NULL) {	/* open input file */
		fprintf (stderr," error: cannot open file '%s'\n", IO_file);
		fprintf (stderr," usage: frame infile\n");
		exit(1);
	}


	lc = 0; /* needs initialising, according to GCC */

	parse_input(fp);
	fclose(fp);

	/*open clean input file*/
	if ((fp = fopen ("frame3dd.cln", "r")) == NULL) {
		fprintf (stderr," error: cannot open cleaned input file \n");
		exit(1);
	}

	frame3dd_getline(fp, title, 256);
	fprintf(stderr," ** %s ** \n\n", title );

	fscanf(fp, "%d %d %d", &nJ, &nM, &nL );
	printf(" number of joints ");
	dots(35);
	printf(" nJ = %d\n", nJ);
	printf(" number of members ");
	dots(34);
	printf(" nM = %d\n", nM);
	printf(" number of load cases ");
	dots(31);
	printf(" nL = %d\n", nL);
        if ( nJ > nM + 1) {
            fprintf(stderr,"warning: %d joints and %d members...", nJ, nM );
	    fprintf(stderr," not enough members to connect all joints.\n");
        }
        if ( nL < 1 ) {
            fprintf(stderr,"error: the number of load cases must be at least 1\n");
	    exit(1);
        }
        if ( nL >= _NL_ ) {
            fprintf(stderr,"error: maximum of %d load cases allowed\n", _NL_-1 );
	    exit(1);
        }

	DoF = 6*nJ;		/* total number of degrees of freedom	*/

				/* allocate memory ... */
//	nF  = ivector(1,nL);	/* # loaded joints, each load case */
//	nW  = ivector(1,nL);	/* # uniformly loaded members,  load case */
//	nP  = ivector(1,nL);	/* # point loaded members, each load case */
//	nT  = ivector(1,nL);	/* # members with temp changes, load case */

	pos = (vec3 *)malloc(sizeof(vec3)*(1+nJ));

	r   = dvector(1,nJ);	/* rigid radius around each joint	*/
	L   = dvector(1,nM);	/* length of each element		*/
	Le  = dvector(1,nM);	/* effective length of each element	*/

	J1  = ivector(1,nM);	/* joint #1 of each element		*/
	J2  = ivector(1,nM);	/* joint #2 of each element		*/
	R   = ivector(1,DoF);	/* reaction force at each degree of freedom */

	Ax  = dvector(1,nM);	/* cross section area of each element	*/
	Asy = dvector(1,nM);	/* shear area in local y direction 	*/
	Asz = dvector(1,nM);	/* shear area in local z direction	*/
	J   = dvector(1,nM);	/* torsional moment of inertia 		*/
	Iy  = dvector(1,nM);	/* bending moment of inertia about local y-axis	*/
	Iz  = dvector(1,nM);	/* bending moment of inertia about local z-axis */
	E   = dvector(1,nM);	/* Young's modulus of elasticity	*/
	G   = dvector(1,nM);	/* shear modulus of elasticity		*/
	p   = dvector(1,nM);	/* member rotation angle about local x axis */

	W   =  D3dmatrix(1,nL,1,nM,1,4); /* distributed load on each member */
	P   =  D3dmatrix(1,nL,1,nM,1,5); /* internal point load each member */
	T   =  D3dmatrix(1,nL,1,nM,1,8); /* internal temp change each member */

	Fo_mech  = dmatrix(1,nL,1,DoF);	/* mechanical load vector	*/
	Fo_temp  = dmatrix(1,nL,1,DoF);	/* temperature load vector	*/
	Fo_temp_lc = dvector(1,DoF); /* external load vector, load case lc */
	Fo_mech_lc = dvector(1,DoF); /* external load vector, load case lc */

	Fo       = dmatrix(1,nL,1,DoF);	/* external load vector		*/
	Fo_lc    = dvector(1,DoF);	/* external load vector, load case lc */
	F_lc     = dvector(1,DoF);	/* external load vector with react'ns */

	feF_mech =  D3dmatrix(1,nL,1,nM,1,12);	/* fixed end forces due to mechanical loads */
	feF_temp =  D3dmatrix(1,nL,1,nM,1,12);	/* fixed end forces due to temperature loads */
	feF      = dmatrix(1,nM,1,12);	/* fixed end forces	*/

	K   = dmatrix(1,DoF,1,DoF);	/* global stiffness matrix	*/
	Q   = dmatrix(1,nM,1,12);	/* end forces for each member	*/

	D   = dvector(1,DoF);	/* displacments of each joint		*/
	dD  = dvector(1,DoF);	/* incremental displacement  of each joint */
	Dp  = dvector(1,DoF);	/* prescribed displacement of each joint */

	d   = dvector(1,nM);	/* mass density for each member		*/
	BMs = dvector(1,nM);	/* lumped beam mass for each member	*/
	JMs = dvector(1,nJ);	/* joint mass for each joint		*/
	JMx = dvector(1,nJ);	/* joint inertia about global X axis	*/
	JMy = dvector(1,nJ);	/* joint inertia about global Y axis	*/
	JMz = dvector(1,nJ);	/* joint inertia about global Z axis	*/

	q = ivector(1,DoF); 	/* vector of condensed degrees of freedom */
	m = ivector(1,DoF); 	/* vector of condensed mode numbers	*/

	read_input_data(
		fp, nJ, nM, pos, r, L, Le, J1, J2, &anlyz, &geom, Q,
		Ax,Asy,Asz, J,Iy,Iz, E,G, p, &shear, mesh_file,plot_file,&exagg
	);
	printf("   input data complete\n");

	assemble_loads (
		fp, nL, nJ, pos, L, Le, Ax,Asy,Asz, Iy,Iz, E, G, p, shear,
		J1, J2, DoF, nM, nF, nW, nP, nT,
		Fo_mech, Fo_temp, W, P, T, feF_mech, feF_temp
	);
	for (i=1; i<=DoF; i++)
		for (lc=1; lc<=nL; lc++)
			Fo[lc][i] = Fo_temp[lc][i] + Fo_mech[lc][i];
	printf("   load data complete\n");

	read_reaction_data ( fp, DoF, &nD, &nR, nJ, Dp, R, &sumR );
	printf("   reaction data complete\n");

	read_mass_data (
		fp, nJ, nM, &nI, d, BMs, JMs, JMx, JMy, JMz, L, Ax,
		&total_mass, &struct_mass, &modes, &Mmethod,
		&lump, mode_file, &tol, &shift, anim, &pan
	);
	printf("   mass data complete\n");

	read_condense ( fp, nJ, modes, &nC, &Cdof, &Cmethod, q, m );

	printf("   matrix condensation data complete\n");

	fclose(fp);
	fp = fopen(IO_file, "a");     /* output appends input */
	if(fp==NULL){
		fprintf(stderr,"Unable to append to input file '%s'!\n",IO_file);
		exit(1);
	}

	write_input_data ( fp, title, nJ,nM,nL, nD,nR, nF,nW,nP,nT,
			pos, r, J1,J2, Ax,Asy,Asz, J,Iy,Iz, E,G, p,
			Fo,Dp,R,W,P,T, shear, anlyz, geom );

	if (anlyz) {			/* solve the problem  */

	 for (lc=1; lc<=nL; lc++) {		/* begin load case loop */

	  for (i=1; i<=DoF; i++)	D[i] = dD[i] = 0.0;
	  for (i=1; i<=DoF; i++)	Fo_temp_lc[i] = Fo_temp[lc][i];
	  for (i=1; i<=DoF; i++)	Fo_mech_lc[i] = Fo_mech[lc][i];

	  assemble_K ( K, DoF, nM, pos,r, L, Le, J1, J2,
			Ax, Asy, Asz, J,Iy,Iz, E, G, p, shear, geom, Q
	  );

#ifdef MATRIX_DEBUG
	  save_dmatrix ( DoF, DoF, K, "Kf" );	     /* free stiffness matrix */
#endif

	  /* apply temperature loads first ... */
	  if (nT[lc] > 0) {
		fprintf(stderr,"\n Linear Elastic Analysis ... Temperature Loads\n");
		apply_reactions ( DoF, R, Dp, Fo_temp_lc, F_lc, K );
		solve_system( K, dD, F_lc, DoF, &ok );
		for (i=1; i<=DoF; i++)	D[i] += dD[i];
		end_forces ( Q, nM, pos, L, Le, J1,J2,
			Ax, Asy,Asz, J,Iy,Iz, E,G, p, D, shear, geom );
	  }

	  assemble_K ( K, DoF, nM, pos, r, L, Le, J1, J2,
			Ax, Asy, Asz, J,Iy,Iz, E, G, p, shear, geom, Q );

	  /* then add mechanical loads ... */
	  if ( nF[lc] > 0 || nW[lc] > 0 || nP[lc] > 0 ) {
		fprintf(stderr,"\n Linear Elastic Analysis ... Mechanical Loads\n");
		apply_reactions ( DoF, R, Dp, Fo_mech_lc, F_lc, K );
		solve_system( K, dD, F_lc, DoF, &ok );
		for (i=1; i<=DoF; i++)	D[i] += dD[i];
		end_forces ( Q, nM, pos, L, Le, J1,J2,
			Ax, Asy,Asz, J,Iy,Iz, E,G, p, D, shear, geom );
	  }

//#ifdef MATRIX_DEBUG
	  save_ut_dmatrix ( DoF, K, "Ks" );	   /* static stiffness matrix */
//#endif


	  for (i=1; i<=DoF; i++)   Fo_lc[i] = Fo_temp_lc[i] + Fo_mech_lc[i];

	  /* Newton-Raphson iterations for geometric non-linearity */
	  if ( geom ) {
		Fe  = dvector( 1, DoF );	/* force equilibrium error  */
		Ks  = dmatrix( 1, DoF, 1, DoF );
		for (i=1;i<=DoF;i++) /* initialize Broyden secant stiffness */
			for(j=i;j<=DoF;j++)
				Ks[i][j]=Ks[j][i]=K[i][j];
	        fprintf(stderr,"\n Non-Linear Elastic Analysis ...\n");
	  }
	  while ( geom && error > tol && iter < 10 && ok >= 0) {

		++iter;

		assemble_K ( K, DoF, nM, pos, r, L, Le, J1, J2,
			Ax, Asy, Asz, J,Iy,Iz, E, G, p, shear, geom, Q );

		apply_reactions ( DoF, R, Dp, Fo_lc, F_lc, K );

		for (i=1; i<=DoF; i++) {	/* equilibrium error	*/
			Fe[i] = F_lc[i];
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
		apply_reactions ( DoF, R, Dp, Fe, Fe, Ks );
		solve_system ( K, dD, Fe, DoF, &ok );

		if ( ok < 0 ) {
		 fprintf(stderr,"   The stiffness matrix is not pos-def. \n");
		 fprintf(stderr,"   Reduce loads and re-run the analysis.\n");
		 break;
		}

		for (i=1; i<=DoF; i++)	D[i] += dD[i];	/* increment D */

		end_forces ( Q, nM, pos, L, Le,
			J1,J2, Ax, Asy,Asz, J,Iy,Iz, E,G, p, D, shear, geom );

						/* convergence criteria: */
#if 0
		error = rel_norm (dD, D, DoF ); /* displacement increment */
#else
		error = rel_norm ( Fe, F_lc, DoF );	/* force balance */
#endif

		fprintf(stderr,"   NR iteration %3d ---", iter);
	        fprintf(stderr," RMS equilibrium precision: %8.2e \n", error);
	  }
	  if ( geom ) {
		free_dvector(Fe, 1, DoF );
		free_dmatrix(Ks, 1, DoF, 1, DoF );
	  }

	  for (i=1; i<=12; i++)
		for (n=1; n<=nM; n++)
			feF[n][i] = feF_temp[lc][n][i] + feF_mech[lc][n][i];
	  equilibrium(pos, L, J1,J2, Fo_lc, R, p, Q, feF, nM, DoF, &error );

	  write_static_results ( fp, nJ,nM,nL,lc, DoF, J1,J2, Fo_lc,
							 D,R,Q, error, ok );

	  write_static_mfile ( argv, nJ,nM,nL,lc, DoF, J1,J2, Fo_lc,
							 D,R,Q, error, ok );

	  static_mesh ( IO_file, mesh_file, plot_file, title, nJ, nM, nL, lc,
				DoF, pos, L, J1,J2, p, D, exagg, anlyz);

	 }					/* end load case loop	*/

	} else {
	    fprintf(stderr,"  %s\n", title );
	    fprintf(stderr,"  data check only\n");
	    static_mesh ( IO_file, mesh_file, plot_file, title, nJ, nM, nL, lc,
				DoF, pos, L, J1,J2, p, D, exagg, anlyz);
	}

	if ( modes > 0 ) {				/* modal analysis */

	        fprintf(stderr,"\n Modal Analysis ...\n");

		calc_modes = (modes+8)<(2*modes) ? modes+8 : 2*modes;

		M   = dmatrix(1,DoF,1,DoF);
		f   = dvector(1,calc_modes);
		V   = dmatrix(1,DoF,1,calc_modes);

		assemble_M ( M, DoF, nJ, nM, pos, r, L, J1, J2, Ax, J, Iy, Iz, p,
					d, BMs, JMs, JMx, JMy, JMz, lump );

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

		save_ut_dmatrix ( DoF, K, "Kd" );	/* dynamic stff matx */
		save_ut_dmatrix ( DoF, M, "Md" );	/* dynamic mass matx */

		if (anlyz) {
		  if ( Mmethod == 1 )
		    subspace( K, M, DoF, calc_modes, f, V, tol,shift,&iter,&ok);
		  if ( Mmethod == 2 )
		    stodola ( K, M, DoF, calc_modes, f, V, tol,shift,&iter,&ok);

		  for (j=1; j<=calc_modes; j++)	f[j] = sqrt(f[j])/(2.*PI);

		  write_modal_results ( fp, nJ,nM,nI, DoF, M,f,V,
				total_mass, struct_mass,
				iter, sumR, modes, shift, lump, tol, ok );
		}
	}

	fclose (fp);

	if ( modes > 0 && anlyz ) {
		modal_mesh ( IO_file, mesh_file, mode_file, plot_file, title,
		       nJ,nM, DoF, modes, pos, L, J1,J2, p, M,f,V,exagg,anlyz);
		animate ( IO_file, mesh_file, mode_file, plot_file, title,anim,
		       nJ,nM, DoF, modes, pos, L, p, J1,J2, f,V, exagg, pan );
	}

	if ( nC > 0 ) {		/* matrix condensation of stiffness and mass */

	        fprintf(stderr,"\n Matrix Condensation ...\n");

		if ( Cdof > modes && Cmethod == 3 ) {
		 fprintf(stderr,"  Cdof > modes ... Cdof = %d  modes = %d \n",
				 Cdof, modes );
		 fprintf(stderr,"  The number of condensed degrees of freedom");
		 fprintf(stderr," may not exceed the number of computed modes");
		 fprintf(stderr," when using dynamic condensation.\n");
		 exit(1);
		}

		Kc = dmatrix(1,Cdof,1,Cdof);
		Mc = dmatrix(1,Cdof,1,Cdof);

		if ( m[1] > 0 && modes > 0 )	Cfreq = f[m[1]];

		if ( Cmethod == 1 ) {	/* static condensation only	*/
			condense ( K, DoF, q, Cdof, Kc);
			printf("   static condensation of K complete\n");
		}
		if ( Cmethod == 2 ) {
			guyan ( M, K, DoF, q, Cdof, Mc,Kc, Cfreq );
			printf("   Guyan condensation of K and M complete");
			printf(" ... dynamics matched at %f Hz.\n", Cfreq );
		}
		if ( Cmethod == 3 && modes > 0 ) {
			dyn_conden ( M,K, DoF, R, q, Cdof, Mc,Kc, V,f, m );
			printf("   dynamic condensation of K and M complete\n");
		}
		save_dmatrix ( Cdof, Cdof, Kc, "Kc" );
		save_dmatrix ( Cdof, Cdof, Mc, "Mc" );

		free_dmatrix ( Kc, 1,Cdof,1,Cdof );
		free_dmatrix ( Mc, 1,Cdof,1,Cdof );
	}

	printf("\n");

	deallocate ( nJ, nM, nL, nF, nW, nP, nT, DoF, modes,
			pos, r, L, Le, J1, J2, R,
			Ax, Asy, Asz, J, Iy, Iz, E, G, p,
			W,P,T,  Fo_mech, Fo_temp, Fo_temp_lc, Fo_mech_lc,
			feF_mech, feF_temp, feF, Fo, Fo_lc, F_lc,
			K, Q, D, dD, Dp,
			d,BMs,JMs,JMx,JMy,JMz, M,f,V, q, m );

	return(0);
}


