/*******************************************************************************

 Program frame.c 
 
 Static and dynamic structural analysis of 2D and 3D frames and trusses with
 elastic and geometric stiffness.
 ---------------------------------------------------------------------------
 http://www.duke.edu/~hpgavin/frame/
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
 Henri P. Gavin                                             hpgavin@duke.edu   
 Department of Civil and Environmental Engineering
 Duke University, Box 90287
 Durham, NC  27708--0287
 ---------------------------------------------------------------------------
 version:    1 March 2007
 ---------------------------------------------------------------------------


Input file format:


A one line descriptive title of your project 

nJ	nM	 				  (number of joints and members)
  J[1]	x[1]	y[1]	z[1]	r[1]	
    :     :       :       :       :		 (joint numbers and coordinates)
  J[nJ]	x[nJ]	y[nJ]	z[nJ]	r[nJ]
					    (member numbers, loc'ns, and prop's)
M[1]  J1[1]  J2[1]  Ax[1]  Asy[1]  Asz[1]  Jp[1]  Iy[1]  Iz[1]  E[1]  G[1]  p[1]
  :      :      :      :       :       :      :      :      :     :     :     ;
M[nM] J1[nM] J2[nM] Ax[nM] Asy[nM] Asz[nM] Jp[nM] Iy[nM] Iz[nM] E[nM] G[nM] p[nM]

shear						  (1: include shear deformation)
geom                              (1: consider geometric nonlinearity, 0: don't)
/tmp/mesh_file						   (mesh data file name)
plot_file					           (mesh plot file name)
exagg						  (exaggerate mesh deformations)
anlyz				     (1: stiffness analysis, 0: data check only)

nF						       (number of loaded joints)
  J[1]	Fx[1]	Fy[1]	Fz[1]	Mxx[1]	Myy[1]	Mzz[1]
   :      :       :       :        :       :       :		  (nodal forces)
  J[nF]	Fx[nF]	Fy[nF]	Fz[nF]	Mxx[nF]	Myy[nF]	Mzz[nF]

nW					   (number of uniform distributed loads)
  M[1]	Wx[1]	Wy[1]	Wz[1]
   :      :       :       :	    (uniform member loads in member coordinates)
  M[nW]	Wx[nW]	Wy[nW]	Wz[nW]

nP					    (number of concentrated point loads)
  M[1]	Px[1]	Py[1]	Pz[1]	x[1]	    (point loads in member coordinates )
    :      :       :       :      :	    (and x=distance from coordinate J1 )
  M[nP]	Px[nP]	Py[nP]	Pz[nP]	x[nP]

nT					 (number of members temperature changes)
  M[1]	a[1]	hy[1]	hz[1]	Ty+[1]	Ty-[1]	Tz+[1]	Tz-[1]	(member no.,   )
    :     :        :       :        :       :       :       :   (temp. coef.)
  M[nT]	a[nT]	hy[nT]	hz[nT]	Ty+[nT]	Ty-[nT]	Tz+[nT]	Tz-[nT] (sizes, & temps)

nR					       (number of joints with reactions)
  J[1]	Rx[1]	Ry[1]	Rz[1]	Rxx[1]	Ryy[1]	Rzz[1]
    :      :       :       :        :       :       :	       (0:free, 1:fixed)
  J[nR]	Rx[nR]	Ry[nR]	Rz[nR]	Rxx[nR]	Ryy[nR]	Rzz[nR]

nD		         (number of joints with prescribed displacements nD<=nR)
  J[1]	Dx[1]	Dy[1]	Dz[1]	Dxx[1]	Dyy[1]	Dzz[1]
    :      :       :       :        :       :       : (prescribed displacements)
  J[nD]	Dx[nD]	Dy[nD]	Dz[nD]	Dxx[nD]	Dyy[nD]	Dzz[nD]

modes						       (number of desired modes)
lump						 (0: consistent mass, 1: lumped)
/tmp/mode_file					     (mode shape data file name)
tol						  (convergence tolerance ~ 1e-4)
shift                                             (eigenvalue shift)

  M[1]	d[1]	BMs[1]
    :     :         :    (beam density and extra beam masses, without self mass)
  M[nM]	d[nM]	BMs[nM]

nI                           (number of joints with extra joint mass or inertia)
  J[1]	JMs[1]    JMx[1]   JMy[1]   JMz[1]  (joint masses and rotatory inertias)
    :       :         :        :        :		    (global coordinates)
  J[nI]	JMs[nI]   JMx[nI]  JMy[nI]  JMz[nI]

nA                                     (number of modes to be animated, nA < 20)
  anim[0] ... anim[nA](list of modes to be animated, sorted by increasing freq.)
pan                                         (1: pan during animation; 0: don't )

nC                                                ( number of condensed joints )
  J[1]  cx[1]  cy[1]  cz[1]   cxx[1]  cyy[1]  czz[1]
    :      :      :      :        :       :       :    ( 1: condense; 0: don't )
  J[nC] cx[nC] cy[nC] cz[nC]  cxx[nC] cyy[nC]  czz[nC]

 ------------------------------------------------------------------------------

to compile:	gcc -O -o frame frame.c eig.c ldl_dcmp.c lu_dcmp.c nrutil.c -lm
to run:		frame 'file-name'

*******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#ifndef PI
#define PI	3.141592653589793
#endif

#define Zvert 1		/* Zvert=0: Y axis is vertical ... rot Z, then rot Y */
			/* Zvert=1: Z axis is vertical ... rot Y, then rot Z */

/* ----------------------- double precision --------------------------------- */
#define float double
#define vector dvector
#define matrix dmatrix
#define free_vector free_dvector
#define free_matrix free_dmatrix
#define save_matrix save_dmatrix
#define show_vector show_dvector 
#define save_ut_matrix save_ut_dmatrix
/* ----------------------- double precision --------------------------------- */
/* also modify all fscanf lines ...     :1,$ s/%f/%lf/g    :1,$ s/%lf/%f/g    */

int main ( argc, argv )
int	argc;
char	*argv[];
{
	char	IO_file[96],	/* the input/output filename		*/
		title[256],	/* the title of the analysis		*/
		mesh_file[96],	/* frame mesh data filename		*/
		plot_file[96],	/* frame mesh plot filename		*/
		mode_file[96],	/* mode-shape mesh data filename	*/
		*strcpy();	/* copy character strings		*/

	FILE	*fp;		/* input/output file pointer		*/

	float	*x, *y, *z,	/* joint coordinates (global)		*/
		*r,		/* joint size radius, for finite sizes	*/
		**K, **Ks,	/* global stiffness matrix		*/
		trK = 0.0,	/* trace of the global stiffness matrix	*/
		**M,		/* global mass matrix			*/
		trM = 0.0,	/* trace of the global mass matrix	*/
		*Fo_mech,	/* mechanical load vector		*/
		*Fo_temp,	/* thermal load vector		*/
		*Fo, *F,	/* a general load vector		*/
		**feF_mech,	/* fixed end forces from mech loads	*/
		**feF_temp,	/* fixed end forces from temp loads	*/
		**feF,		/* a general set of fixed end forces	*/
		*D, *dD,	/* displacement and displ increment	*/
		dDdD = 0.0,	/* dD' * dD				*/
		*Fe,		/* equilibrium error in nonlinear anlys	*/
		*Dp,		/* prescribed joint displacements	*/
		**W,		/* uniform distributed member loads	*/
		**P,		/* member concentrated loads		*/ 
		**T,		/* member temperature  loads		*/
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
		shift,		/* shift-factor for rigid-body-modes	*/
		*f, **V,	/* natural frequencies and mode-shapes	*/
		error = 1.0,	/* rms equilibrium error and reactions	*/
		**Kc, **Mc,	/* condensed stiffness and mass matrices*/
		exagg,		/* exaggerate deformations in mesh data	*/
		rel_norm(),	/* relative 2-norm between two  vectors	*/
		*vector(),	/* dynamic memory allocation for vector	*/
		**matrix();	/* dynamic memory allocation for matrix */

	int	nJ, nM,		/* number of Joints and Members		*/
		DoF, i, j, n,	/* number of Degrees of Freedom		*/
		nF, nR,		/* number of loaded and restrained joints*/
		nD,		/* number of prescribed nodal displ'nts	*/
		nW,		/* number of members w/ unif distr loads*/ 
		nP,		/* number of members w/ conc point loads*/
		nT,		/* number of members w/ temp. changes	*/
		nI,             /* number of joints w/ extra inertia	*/
		nC,		/* number of condensed joints		*/
		*J1, *J2,	/* begin and end joint numbers		*/
		shear,		/* indicates shear deformationi		*/
		geom=0,		/* indicates  geometric nonlinearity	*/
		anlyz=1,	/* 1: stiffness analysis, 0: data check	*/
		*R, sumR,	/* reaction data, total no. of reactions*/
		modes,		/* number of desired modes		*/
		calc_modes,	/* number of modes to calculate		*/
		lump=1,		/* 1: lumped, 0: consistent mass matrix */
		iter=0,		/* number of iterations			*/
		ok=1,		/* number of (-ve) diag. terms of L D L' */
		anim[20],	/* the modes to be animated		*/
		pan=1,		/* 1: pan during animation; 0: don't	*/         
		*q,		/* vector of DoF's to condense		*/
		Cdof,		/* number of condensed degrees o freedom*/
		temp_mech,	/* counter for temp and mech load cases	*/
		*ivector();

	void	parse_input(),	/* re-write input file without comments */
		read_input(),	/* read input data file			*/
		getline(),	/* get a line from the input file	*/
		assemble_K(),	/* form the global stiffness matrix	*/
		read_loads(),	/* form load vector			*/
		read_reactions(),  /* read boundary conditions		*/
		apply_reactions(), /* apply boundary conditions		*/
		solve_system(),	/* solve a linear system via LDL' dcmp	*/
		end_forces(),	/* evaluate the member end forces	*/
		equilibrium(),	/* perform an equilibrium check		*/
		read_masses(),	/* read density and extra inertial mass	*/
		assemble_M(),	/* form the global mass matrix		*/
		read_condense(),/* read matrix condensation data	*/
		condense(),	/* static matrix condensation		*/
		guyan(),	/* Guyan reduction of matrices Md , Kd	*/
		dyn_conden(),	/* dynamic condensation of Md and Kd	*/
		save_matrix(),	/* save a matrix (for debugging)	*/
		save_ut_matrix(),/* save a matrix (for debugging)	*/
		save_ivector(),	/* save a vector of integers		*/
		stodola(),	/* lower generalized eigenval & eigenvec*/
		subspace(),	/* lower generalized eigenval & eigenvec*/
		control_data(),	/* save control data information	*/
		save_results(),	/* save displacements and member forces	*/
		mesh(),		/* create undeformed and deformed meshes*/
		modal_results(),/* save nat'l frequencies & mode shapes	*/
		modal_mesh(),	/* create undeformed and mode-shape meshes*/
		animate(),	/* animate the mode shape meshes	*/
		dots(),		/* print a string of dots (periods)	*/
		free_vector(),	/* deallocate memory in a vector	*/
		free_matrix(),	/* deallocate memory in a matrix	*/
		deallocate(),	/* release the allocated memory		*/
		exit();		/* exit the program			*/


        fprintf(stderr," FRAME version:  1 Mar 2007,");
        fprintf(stderr," GPL Copyright (C) 1992-2007, Henri P. Gavin \n");
        fprintf(stderr," http://www.duke.edu/~hpgavin/frame/ \n");
	fprintf(stderr," This is free software with absolutely no warranty.\n");
	fprintf(stderr," For details, see http://www.fsf.org/copyleft/gpl.html \n\n");

	if (argc < 2) {
		fprintf (stderr," Please enter the input/output file name: ");
		scanf("%s", IO_file );
		fprintf (stderr," You entered file name: %s \n", IO_file );
	} else strcpy ( IO_file , argv[1] );

	if ((fp = fopen (IO_file, "r")) == NULL) {	/* open input file */
		fprintf (stderr," error: cannot open file '%s'\n", IO_file);
		fprintf (stderr," usage: frame infile\n");
		exit(1);
	}

	parse_input(fp);
	fclose(fp);

	if ((fp = fopen ("frame.cln", "r")) == NULL) { /*open clean input file*/
		fprintf (stderr," error: cannot open cleaned input file \n");
		exit(1);
	}

	getline(fp, title, 256);
	fprintf(stderr," ** %s ** \n\n", title );

	fscanf(fp, "%d %d", &nJ, &nM );
	printf(" number of joints ");
	dots(35);
	printf(" nJ = %d\n", nJ);
	printf(" number of members ");
	dots(34);
	printf(" nM = %d\n", nM);
        if ( nJ > nM + 1) {
                fprintf(stderr,"warning: %d joints and %d members...", nJ, nM );
		fprintf(stderr," not enough members to connect all joints.\n");
        }

	DoF = 6*nJ;		/* total number of degrees of freedom		*/

							   /* allocate memory */
	x   =  vector(1,nJ);	/* x coordinate of each joint		*/
	y   =  vector(1,nJ);	/* y coordinate of each joint		*/
	z   =  vector(1,nJ);	/* z coordinate of each joint		*/
	r   =  vector(1,nJ);	/* rigid radius around each joint	*/
	L   =  vector(1,nM);	/* length of each element		*/
	Le  =  vector(1,nM);	/* effective length of each element	*/
	J1  = ivector(1,nM);	/* joint #1 of each element		*/
	J2  = ivector(1,nM);	/* joint #2 of each element		*/
	Ax  =  vector(1,nM);	/* cross section area of each element	*/
	Asy =  vector(1,nM);	/* shear area in local y direction 	*/
	Asz =  vector(1,nM);	/* shear area in local z direction	*/
	J   =  vector(1,nM);	/* torsional moment of inertia 		*/
	Iy  =  vector(1,nM);	/* bending moment of inertia about local y-axis	*/
	Iz  =  vector(1,nM);	/* bending moment of inertia about local z-axis */
	E   =  vector(1,nM);	/* Young's modulus of elasticity	*/
	G   =  vector(1,nM);	/* shear modulus of elasticity		*/
	p   =  vector(1,nM);	/* member rotation angle about local x axis */
	d   =  vector(1,nM);	/* mass density for each member		*/
	BMs =  vector(1,nM);	/* lumped beam mass for each member	*/
	JMs =  vector(1,nJ);	/* joint mass for each joint		*/
	JMx =  vector(1,nJ);	/* joint inertia about global X axis	*/
	JMy =  vector(1,nJ);	/* joint inertia about global Y axis	*/
	JMz =  vector(1,nJ);	/* joint inertia about global Z axis	*/

	K   =  matrix(1,DoF,1,DoF);	/* global stiffness matrix	*/
	Q   =  matrix(1,nM,1,12);	/* end forces for each member	*/

	D   =  vector(1,DoF);	/* displacments of each joint		*/
	dD  =  vector(1,DoF);	/* incremental displacement  of each joint */
	Dp  =  vector(1,DoF);	/* prescribed displacement of each joint */
	R   = ivector(1,DoF);	/* reaction force at each degree of freedom */
	W   =  matrix(1,nM,1,4);	/* distributed load on each member */
	P   =  matrix(1,nM,1,5);	/* internal point load each member */
	T   =  matrix(1,nM,1,8);	/* internal temp change each member */
	Fo_mech  =  vector(1,DoF);	/* mechanical load vector	*/
	Fo_temp  =  vector(1,DoF);	/* temperature load vector	*/
	Fo       =  vector(1,DoF);	/* external load vector		*/
	F        =  vector(1,DoF);	/* external load vector inc'dng reactions */
	feF_mech =  matrix(1,nM,1,12);	/* fixed end forces due to mechanical loads */
	feF_temp =  matrix(1,nM,1,12);	/* fixed end forces due to temperature loads */
	feF      =  matrix(1,nM,1,12);	/* fixed end forces	*/

	q = ivector(1,DoF); 	/* vector of condensed degrees of freedom */

	read_input ( fp, nJ, nM, x,y,z,r, L, Le, J1, J2, &anlyz, &geom, Q, 
	     Ax,Asy,Asz, J,Iy,Iz, E,G, p, &shear, mesh_file,plot_file,&exagg);
	printf("read input done\n");

	read_loads ( fp, nJ, x, y, z, L, Le, Ax,Asy,Asz, Iy,Iz, E, G, p, shear, 
		 J1, J2, DoF, nM, &nF, &nW, &nP, &nT,
		 Fo_mech, Fo_temp, W, P, T, feF_mech, feF_temp );
	for (i=1; i<=DoF; i++)	Fo[i] = Fo_temp[i] + Fo_mech[i];
	printf("read loads done\n");

	read_reactions ( fp, DoF, &nD, &nR, nJ, Dp, R, &sumR );
	printf("read reactions done\n");

	read_masses ( fp, nJ, nM, &nI, d, BMs, JMs, JMx, JMy, JMz, L, Ax, 
			&total_mass, &struct_mass, &modes, &lump, mode_file,
			&tol, &shift, anim, &pan );
	printf("read masses done\n");

	read_condense ( fp, nJ, &nC, &Cdof, q );

	fclose (fp);	fp = fopen(IO_file, "a");     /* output appends input */

	if(fp==NULL){
		fprintf(stderr,"Unable to open input file for writing (appending output data)\n");
		exit(1);
	}
 
	control_data ( fp, title, nJ,nM, nF,nD,nR,nW,nP,nT, x,y,z,r, J1,J2, 
		Ax,Asy,Asz, J,Iy,Iz, E,G, p, Fo,Dp,R,W,P,T, shear,anlyz,geom );

	if (anlyz) {					/* solve the problem  */

	  for (i=1; i<=DoF; i++)	D[i] = dD[i] = 0.0;

	  for (temp_mech=1; temp_mech <= 2; temp_mech++) { 

	    assemble_K ( K, DoF, nM, x,y,z,r, L, Le, J1, J2,
			Ax, Asy, Asz, J,Iy,Iz, E, G, p, shear, geom, Q );

/*	    save_matrix ( DoF, DoF, K, "Kf" );	     /* free stiffness matrix */

	    if (temp_mech == 1 && nT > 0) { /* temperature loads first	*/
		fprintf(stderr,"\n Linear Elastic Analysis ... Temperature Loads\n");
		apply_reactions ( DoF, R, Dp, Fo_temp, F, K );
		solve_system( K, dD, F, DoF, &ok );
		for (i=1; i<=DoF; i++)	D[i] += dD[i];
		end_forces ( Q, nM, x,y,z, L, Le, J1,J2, Ax, Asy,Asz, J,Iy,Iz, E,G, p, D, shear,geom );
	    }
	    if (temp_mech == 2 && (nF > 0 || nW > 0 || nP > 0) ) { /* then add mechanical loads	*/
		fprintf(stderr,"\n Linear Elastic Analysis ... Mechanical Loads\n");
		apply_reactions ( DoF, R, Dp, Fo_mech, F, K );
		solve_system( K, dD, F, DoF, &ok );
		for (i=1; i<=DoF; i++)	D[i] += dD[i];
		end_forces ( Q, nM, x,y,z, L, Le, J1,J2, Ax, Asy,Asz, J,Iy,Iz, E,G, p, D, shear,geom );
	    } 

/*	    save_matrix ( DoF, DoF, K, "Ks" );	   /* static stiffness matrix */

	  } 

	  for (i=1; i<=DoF; i++)	Fo[i] = Fo_temp[i] + Fo_mech[i];

	  /* Newton-Raphson iterations for geometric non-linearity */
	  if ( geom ) {
		Fe  =  vector( 1, DoF );	/* force equilibrium error  */
		Ks  =  matrix( 1, DoF, 1, DoF ); 
		for (i=1;i<=DoF;i++) /* initialize Broyden secant stiffness */
			for(j=i;j<=DoF;j++)
				Ks[i][j]=Ks[j][i]=K[i][j];
	        fprintf(stderr," Non-Linear Elastic Analysis ...\n");
	  }
	  while ( geom && error > tol && iter < 10 && ok >= 0) { 

		++iter;

		assemble_K ( K, DoF, nM, x,y,z,r, L, Le, J1, J2,
			Ax, Asy, Asz, J,Iy,Iz, E, G, p, shear, geom, Q );

		apply_reactions ( DoF, R, Dp, Fo, F, K );

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
		apply_reactions ( DoF, R, Dp, Fe, Fe, Ks );
		solve_system ( K, dD, Fe, DoF, &ok );
		 
		if ( ok < 0 ) {
		 fprintf(stderr,"   The stiffness matrix is not pos-def. \n");
		 fprintf(stderr,"   Reduce loads and re-run the analysis.\n");
		 break;
		}

		for (i=1; i<=DoF; i++)	D[i] += dD[i];	/* increment D */

		end_forces ( Q, nM, x,y,z, L, Le,
			J1,J2, Ax, Asy,Asz, J,Iy,Iz, E,G, p, D, shear, geom );

						/* convergence criteria: */
		/* error = rel_norm (dD, D, DoF ); /* displacement increment */
		error = rel_norm ( Fe, F, DoF );	/* force balance */

		fprintf(stderr,"  ... NR iter %2d:  error = %8.2e \n",
								iter, error );
	  }
	  if ( geom ) {
		free_vector(dD, 1, DoF );
		free_vector(Fe, 1, DoF );
		free_matrix(Ks, 1, DoF, 1, DoF );
	  }
	    
	  for (i=1; i<=12; i++)
		for (n=1; n<=nM; n++)
			feF[n][i] = feF_temp[n][i] + feF_mech[n][i];
	  equilibrium ( x,y,z, L, J1,J2, Fo, R, p, Q, feF, nM, DoF, &error );
	  save_results ( fp, nJ, nM, DoF, J1, J2, Fo, D, R, Q, error, ok );

	} else {
	    fprintf(stderr,"  %s\n", title );
	    fprintf(stderr,"  data check only\n");
	}

	if ( modes > 0 ) {

		calc_modes = (modes+8)<(2*modes) ? modes+8 : 2*modes;

		M   =  matrix(1,DoF,1,DoF);
		f   =  vector(1,calc_modes);
		V   =  matrix(1,DoF,1,calc_modes);

		assemble_M ( M, DoF, nJ, nM, x,y,z,r, L, J1, J2, Ax, J, Iy, Iz, p,
					d, BMs, JMs, JMx, JMy, JMz, lump );

/*		save_matrix ( DoF, DoF, M, "Mf" );	/* free mass matrix */

		for (j=1; j<=DoF; j++) { trK += K[j][j]; trM += M[j][j]; }
		for (i=1; i<=DoF; i++) {
			if ( R[i] ) {	/* apply reactions to upper triangle */
				K[i][i] = trK * 1e2; 	
				M[i][i] = trM; 
				for (j=i+1; j<=DoF; j++) 
					K[j][i]=K[i][j]=M[j][i]=M[i][j] = 0.0;
		    }
		}

		save_ut_matrix ( DoF, K, "Kd" );	/* dynamic stff matx */
		save_ut_matrix ( DoF, M, "Md" );	/* dynamic mass matx */
	
		if (anlyz) {
		  if (0)
		    stodola ( K, M, DoF, calc_modes, f, V, tol,shift,&iter,&ok);
		  if (1)
		    subspace( K, M, DoF, calc_modes, f, V, tol,shift,&iter,&ok);
	
		  for (j=1; j<=calc_modes; j++)	f[j] = sqrt(f[j])/(2.*PI);

		  modal_results ( fp, nJ,nM,nI, DoF, M,f,V,
				total_mass, struct_mass,
				iter, sumR, modes, shift, lump, tol, ok );
		}
	}

	fclose (fp);

	mesh ( IO_file, mesh_file, plot_file, title, nJ, nM, DoF, 
					x,y,z,L, J1,J2, p, D, exagg, anlyz);

	if ( modes > 0 && anlyz ) {
		modal_mesh ( IO_file, mesh_file, mode_file, plot_file, title, 
		       nJ,nM, DoF, modes, x,y,z,L, J1,J2, p, M,f,V,exagg,anlyz);
		animate ( IO_file, mesh_file, mode_file, plot_file, title,anim, 
		       nJ,nM, DoF, modes, x,y,z,L, p, J1,J2, f,V, exagg, pan );
	}

	if ( nC > 0 ) {		/* dynamic condensation of stiffness and mass */
		Kc = matrix(1,Cdof,1,Cdof);
		Mc = matrix(1,Cdof,1,Cdof);
		condense ( K, DoF, q, Cdof, Kc);
		save_matrix ( Cdof, Cdof, Kc, "Kcc" );	
		if ( modes > 0 && 1 )
			guyan ( M, K, DoF, q, Cdof, Mc,Kc, f[1] );
		if ( modes > 0 && 0 )
			dyn_conden ( M,K, DoF, q, Cdof, Mc,Kc, V,f );
		save_matrix ( Cdof, Cdof, Kc, "Kc" );	
		save_matrix ( Cdof, Cdof, Mc, "Mc" );  

		free_matrix ( Kc, 1,Cdof,1,Cdof );
		free_matrix ( Mc, 1,Cdof,1,Cdof );
	}

	deallocate ( x,y,z,r, L,Le, J1, J2, Ax,Asy,Asz, J,Iy,Iz, E, G,
       K,Q,F,D,R,W,P,T, feF,Fo, d,BMs,JMs,JMx,JMy,JMz,M,f,V, nJ,nM,DoF, modes );

	return(0);
}


/*------------------------------------------------------------------------------
READ_INPUT  -  read material and geometry data, calc lengths		15dec97
------------------------------------------------------------------------------*/
void read_input ( fp, nJ, nM, x,y,z,r, L,Le, J1, J2, anlyz, geom, Q, 
		Ax,Asy,Asz, J,Iy,Iz, E,G, p, shear, meshfile, plotfile, exagg )
char	meshfile[], plotfile[];
FILE	*fp;
float	*x, *y, *z, *r, *L, *Le, *Ax,*Asy,*Asz, *J, *Iy,*Iz, *E, *G, *p, *exagg;
float	**Q;
int	nJ, nM, *J1, *J2, *shear, *geom, *anlyz; 
{
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
void getline (fp, s, lim)
FILE	*fp;
char    *s;
int     lim;
{
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
void parse_input(fp)
FILE *fp;
{
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
void getline_no_comment (fp, s, lim)
FILE    *fp;            /* pointer to the file from which to read */
char    *s;             /* pointer to the string to which to write */
int     lim;            /* the longest anticipated line length  */
{
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
void read_loads ( fp, nJ, x, y, z, L, Le, Ax, Asy, Asz, Iy, Iz, E, G, p, shear, 
		 J1, J2, DoF, nM, nF, nW, nP, nT, 
		 F_mech, F_temp, W, P, T, feF_mech, feF_temp )
FILE	*fp;
int	nJ, *J1, *J2, DoF, *nF, *nW, *nP, *nT, nM, shear;
float	*x, *y, *z, *L, *Le, *Ax, *Asy, *Asz, *Iy, *Iz, *E, *G, *p, 
	*F_mech, *F_temp, **W, **P, **T, **feF_mech, **feF_temp;
{
	float	Nx1, Vy1, Vz1, Mx1, My1, Mz1,	/* fixed end forces */
		Nx2, Vy2, Vz2, Mx2, My2, Mz2,
		Ksy, Ksz, 		/* shear deformatn coefficients	*/
		a, b,				/* point load locations */
		hy, hz,			/* section dimensions in local coords */
		t1, t2, t3, t4, t5, t6, t7, t8, t9;	/* 3D coord Xfrm coef */
	int	i,j,l,n, j1, j2;
	void	coord_trans(),	/* 3D coordinate transformation	*/
		dots(),
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
void read_reactions ( fp, DoF, nD, nR, nJ, Dp, R, sumR )
FILE	*fp;
int	DoF, *nD, *nR, nJ, *R, *sumR;
float	*Dp;
{
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
ASSEMBLE_K  -  assemble global stiffness matrix from individual elements 23feb94
------------------------------------------------------------------------------*/
void assemble_K ( K, DoF,nM, x,y,z,r, L,Le, J1,J2, Ax,Asy,Asz, J,Iy,Iz, E,G, p,
			shear, geom, Q )
int	DoF, nM;
float	**K, *x,*y,*z,*r, *L,*Le, *Ax, *Asy,*Asz, *J, *Iy,*Iz, *E, *G, *p, **Q;
int	*J1, *J2, shear, geom;
{
	float	**k,		/* element stiffness matrix in global coord */
		**matrix();
	int	**ind,		/* member-structure DoF index table	*/
		**imatrix(),
		i, j, ii, jj, l, ll;
	void	elastic_K(),
		geometric_K(),
		end_release(), 
		save_matrix(), free_matrix(), free_imatrix();

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
		  Asy[i],Asz[i], Iy[i],Iz[i], E[i],G[i], p[i], -Q[i][1], shear);

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
void elastic_K( k, x,y,z,r, L,Le, j1, j2, Ax, Asy,Asz, J, Iy,Iz, E,G, p, shear )
float	**k, *x, *y, *z, *r, L, Le, Ax, Asy, Asz, J, Iy, Iz, E, G, p;
int     j1, j2, shear;
{
	float   t1, t2, t3, t4, t5, t6, t7, t8, t9,     /* coord Xformn */
		Ksy, Ksz;		/* shear deformatn coefficients	*/
	int     i, j;
	void    coord_trans(),  /* 3D coordinate transformation coeff's */
		save_matrix(),
		atma();		/* carry out the coordinate transfm'n	*/

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
GEOMETRIC_K - space frame geometric stiffness matrix, global coordnates 12nov02
------------------------------------------------------------------------------*/
void geometric_K( k, x,y,z,r, L, Le, j1, j2, Asy, Asz, Iy,Iz,E,G, p, T, shear )
float	**k, *x, *y, *z, *r, L, Le, Asy, Asz, Iy, Iz, E, G, p, T;
int     j1, j2, shear;
{
	float   t1, t2, t3, t4, t5, t6, t7, t8, t9,     /* coord Xformn */
		**kg,
		Ksy,Ksz,Dsy,Dsz,/* shear deformation coefficients	*/
		**matrix();
	int     i, j;
	void    coord_trans(),  /* 3D coordinate transformation coeff's */
		save_matrix(),
		free_matrix(),
		atma();		/* carry out the coordinate transfm'n	*/


	coord_trans ( x, y, z, L, j1, j2,
			&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p );

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


	kg[2][2]  = kg[8][8]   = T/L*(1.2+2.0*Ksy+Ksy*Ksy)/Dsy;
	kg[3][3]  = kg[9][9]   = T/L*(1.2+2.0*Ksz+Ksz*Ksz)/Dsz;
	kg[5][5]  = kg[11][11] = T*L*(2.0/15.0+Ksz/6.0+Ksz*Ksz/12.0)/Dsz;
	kg[6][6]  = kg[12][12] = T*L*(2.0/15.0+Ksy/6.0+Ksy*Ksy/12.0)/Dsy;

	kg[5][3]  = kg[3][5]   =  kg[11][3] = kg[3][11] = -T/10.0/Dsz;
	kg[9][5]  = kg[5][9]   =  kg[11][9] = kg[9][11] =  T/10.0/Dsz;
	kg[6][2]  = kg[2][6]   =  kg[12][2] = kg[2][12] =  T/10.0/Dsy;
	kg[8][6]  = kg[6][8]   =  kg[12][8] = kg[8][12] = -T/10.0/Dsy;

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
		error=1.0,	/* error in the solution		*/
		*vector();	/* vector memory allocation		*/

	void	ldl_dcmp(),	/* L D L' decompositon and back-sub'n	*/
		ldl_mprove(),	/* iterative improvement to the solution  */
		free_vector(),	/* vector memory de-allocation		*/
		exit();

	diag = vector ( 1, DoF );

	ldl_dcmp( K, DoF, diag, F, D, 1, 0, ok );	/*  L D L'  decomp */
	if ( *ok < 0 ) {
	 	fprintf(stderr," Make sure that all six");
		fprintf(stderr," rigid body translations are restrained!\n");
		/* exit(1); */
	} else {				/* back substitute for D */
		ldl_dcmp( K, DoF, diag, F, D, 0, 1, ok ); /* LDL'  back-sub */
	        fprintf(stderr,"  RMS matrix error:"); /* improve solution */
		error = *ok = 1;
		do {
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
	float	*s, *vector();
	int	i,j;
	void	member_force(),
		free_vector();

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
void member_force ( s, M, x,y,z, L, Le,
			j1,j2, Ax, Asy,Asz, J, Iy,Iz, E,G,p, D, shear, geom )
int	M, j1, j2, shear, geom;
float	*s, *x, *y, *z, L, Le, Ax, Asy, Asz, J, Iy, Iz, E, G, p, *D; 
{
	float	t1, t2, t3, t4, t5, t6, t7, t8, t9, /* coord Xformn	*/
		d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12,
		x1, y1, z1, x2, y2, z2,	/* joint coordinates	*/
		Ls,			/* stretched length of element */
		Ksy, Ksz, Dsy, Dsz,	/* shear deformation coeff's	*/
		T;		/* normal force for geometric stiffness */
	void	coord_trans(); /* 3D coordinate transformation	*/

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
	void	coord_trans(); /* 3D coordinate transformation coeff's */

	for (j=1; j<=DoF; j++)	if ( R[j] == 0 )   den += ( F[j]*F[j] );
	den = sqrt ( den / (float) DoF );

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
	fprintf(stderr,"  RMS relative equilibrium error: %9.3e\n", *err );

	return;
}


/*------------------------------------------------------------------------------
READ_MASSES  -  read member densities and extra inertial mass data	16aug01
------------------------------------------------------------------------------*/
void read_masses ( fp, nJ, nM, nI, d, BMs, JMs, JMx, JMy, JMz, L, Ax,
	total_mass, struct_mass, modes, lump, modefile, tol,shift,anim, pan)
char	modefile[];
FILE	*fp;
float	*d, *BMs, *JMs, *JMx, *JMy, *JMz, *L, *Ax, *tol, *shift,
	*total_mass, *struct_mass;
int	nJ, nM, *nI, *modes, *lump, anim[], *pan;
{
	float	ms = 0.0;
	int	chk, j, jnt, m, mem, nA;
	void	dots(), exit();

	*total_mass = *struct_mass = 0.0;	

	chk = fscanf ( fp, "%d", modes );
	printf(" number of dynamic modes ");
	dots(28);
	printf(" modes = %d\n",*modes);
	if ( *modes < 1 || chk != 1 ) {
		*modes = 0;
		return;
	} else {
		fscanf( fp, "%d", lump );
		fscanf( fp, "%s", modefile );
		fscanf( fp, "%lf", tol );
		fscanf( fp, "%lf", shift );
		for (m=1; m <= nM; m++) {	/* read inertia data	*/
			fscanf(fp, "%d", &mem );
			fscanf(fp, "%lf %lf", &d[mem], &BMs[mem] );
			*total_mass  += d[mem]*Ax[mem]*L[mem] + BMs[mem];
			*struct_mass += d[mem]*Ax[mem]*L[mem];
		}

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
	void    lumped_M(), consistent_M(),
		save_matrix(), free_matrix(), free_imatrix();

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
void lumped_M ( m, x,y,z, L, j1, j2, Ax, J,Iy,Iz, d, p, BMs )
float   **m, *x, *y, *z, L, Ax, J,Iy,Iz, d, BMs, p;
int     j1, j2;
{
	float   t1, t2, t3, t4, t5, t6, t7, t8, t9,     /* coord Xformn */
		t, ry,rz, po;	/* translational, rotational & polar inertia */
	int     i, j;
	void    coord_trans();  /* 3D coordinate transformation coeff's */

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
void consistent_M ( m, x,y,z,r, L, j1, j2, Ax, J, Iy, Iz, d, BMs, p )
float   **m, *x, *y, *z, *r, L, Ax, J, Iy, Iz, d, BMs, p;
int     j1, j2;
{
	float   t1, t2, t3, t4, t5, t6, t7, t8, t9,     /* coord Xformn */
		t, ry, rz, po;	/* translational, rotational & polar inertia */
	int     i, j;
	void    coord_trans(),  /* 3D coordinate transformation coeff's */
		save_matrix(),
		atma();		/* carry out the coordinate transfm'n	*/

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
COORD_TRANS -  evaluate the 3D coordinate transformation coefficients 1dec04
Default order of coordinate rotations ...  typical for Y as the vertical axis
1. rotate about the global Z axis
2. rotate about the global Y axis
3. rotate about the local  x axis --- element 'roll'

If Zvert is defined as 1, then the order of coordinate rotations is typical
for Z as the vertical axis
1. rotate about the global Y axis
2. rotate about the global Z axis
3. rotate about the local  x axis --- element 'roll'

    Q=TF;   U=TD;   T'T=I;   Q=kU;   TF=kTD;   T'TF=T'kTD;   T'kT = K;   F=KD
------------------------------------------------------------------------------*/
void coord_trans ( x, y, z, L, j1, j2, t1,t2,t3, t4,t5,t6, t7,t8,t9, p )
int	j1, j2;
float	*x, *y, *z, L, 
	*t1, *t2, *t3, *t4, *t5, *t6, *t7, *t8, *t9,
	p;					/* the roll angle (radians) */
{
	float	Cx, Cy, Cz, den,		/* direction cosines	*/
		Cp, Sp;			/* cosine and sine of roll angle */

	Cx = (x[j2] - x[j1]) / L;
	Cy = (y[j2] - y[j1]) / L;
	Cz = (z[j2] - z[j1]) / L;

	*t1 = *t2 = *t3 = *t4 = *t5 = *t6 = *t7 = *t8 = *t9 = 0.0;

	Cp = cos(p);
	Sp = sin(p);

#if Zvert

	if ( fabs(Cz) == 1.0 ) {
		*t3 =  Cz;
		*t4 = -Cz*Sp;
		*t5 =  Cp;
		*t7 = -Cz*Cp;
		*t8 = -Sp;
	} else {

		den = sqrt ( 1.0 - Cz*Cz );

		*t1 = Cx;
	   	*t2 = Cy;
		*t3 = Cz;

		*t4 = (-Cx*Cz*Sp - Cy*Cp)/den;    
 		*t5 = (-Cy*Cz*Sp + Cx*Cp)/den;
		*t6 = Sp*den;

		*t7 = (-Cx*Cz*Cp + Cy*Sp)/den;
		*t8 = (-Cy*Cz*Cp - Cx*Sp)/den;
	   	*t9 = Cp*den;
	}

#else

	if ( fabs(Cy) == 1.0 ) {
		*t2 =  Cy;
		*t4 = -Cy*Cp;
		*t6 =  Sp;
		*t7 =  Cy*Sp;
		*t9 =  Cp;
	} else {

		den = sqrt ( 1.0 - Cy*Cy );

		*t1 = Cx;
	   	*t2 = Cy;
		*t3 = Cz;

		*t4 = (-Cx*Cy*Cp - Cz*Sp)/den;    
		*t5 = den*Cp;
 		*t6 = (-Cy*Cz*Cp + Cx*Sp)/den;

		*t7 = (Cx*Cy*Sp - Cz*Cp)/den;
	   	*t8 = -den*Sp;
		*t9 = (Cy*Cz*Sp + Cx*Cp)/den;
	}

#endif

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
	void	save_matrix(), free_matrix();

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
READ_CONDENSE   -  read matrix condensation information 	        30aug01
------------------------------------------------------------------------------*/
void read_condense ( fp, nJ, nC, Cdof, q )
FILE	*fp;
int	nJ, *nC, *Cdof, *q;
{
	int	i,j,k,  chk, **qm, *ivector(), **imatrix();
	void	dots(),
		free_imatrix(), exit();

	*Cdof = 0;

	if ( (chk = fscanf ( fp, "%d", nC )) != 1 )  {
		*nC = 0;
		return;
	}

	if ( *nC <= 0 ) {
		*nC = 0;
		return;
	}
	printf(" number of joints with condensed DoF's ");
	dots(14);
	printf(" nC = %d\n",*nC);

	qm = imatrix( 1, *nC, 1,7 );

	for ( i=1; i <= *nC; i++) {
	 fscanf( fp, "%d %d %d %d %d %d %d",
	 &qm[i][1],
	 &qm[i][2], &qm[i][3], &qm[i][4], &qm[i][5], &qm[i][6], &qm[i][7]);
	 if ( qm[i][1] < 1 || qm[i][1] > nJ ) {		/* error check */
	  fprintf(stderr," error in matrix condensation data: condensed joint number out of range\n");
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

	free_imatrix(qm,1, *nC, 1,7);
	return;
}


/*------------------------------------------------------------------------------
CONDENSE - static condensation of stiffness matrix from NxN to nxn    30aug01
------------------------------------------------------------------------------*/
void condense ( A, N, q, n, Ac)
float	**A, **Ac;
int	N, n, *q;
{
	float	**Arr, **Arq, **matrix();
	int	i,j,k, ri,rj,qi,qj, ok, 
		*r, *ivector();
	void	xtaiy(),
		save_matrix(), save_ivector(), save_ut_matrix(),
		free_matrix(), free_ivector();


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

	xtaiy ( Arq, Arr, Arq, N-n, n, Ac );

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
GUYAN  -   generalized Guyan reduction of mass and stiffness matrices    7oct01
           matches the response at a particular frequency, sqrt(L)/2/pi
           Guyan, Robert J., ``Reduction of Stiffness and Mass Matrices,''
           AIAA Journal, Vol. 3, No. 2 (1965) p 380.
-----------------------------------------------------------------------------*/
void guyan ( M, K, N, q, n, Mc, Kc, L )
float	**M, **K, **Mc, **Kc, L;
int	N, n, *q;
{
	float	**Drr, **Drq, **Mrr, **Mrq, **Ac, **matrix();
	int	i,j,k, ri,rj,qi,qj, ok, 
		*r, *ivector();
	void	xtaiy(), /* compute X' * inv(A) * Y	*/
		aixai(), /* compute inv(A) * X * inv(A)	*/
		xAx(),	 /* compute X' * A * X	*/
		save_matrix(), save_ivector(), save_ut_matrix(),
		free_matrix(), free_ivector();

	r   = ivector(1,N-n);
	Drr =  matrix(1,N-n,1,N-n);
	Drq =  matrix(1,N-n,1,n);
	Mrr =  matrix(1,N-n,1,N-n);
	Mrq =  matrix(1,N-n,1,n);
	Ac  =  matrix(1,n,1,n);

	L = 4.0 * PI * PI * L * L;	/* eigen-value	*/

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
		for (j=i; j<=N-n; j++) { /* use only upper triangle of K,M */
			ri = r[i];
			rj = r[j];
			if ( ri <= rj )	{
				Mrr[j][i] = Mrr[i][j] = M[ri][rj];
				Drr[j][i] = Drr[i][j] = K[ri][rj] - L*M[ri][rj];
			}
		}
	}

	for (i=1; i<=N-n; i++) {
		for (j=1; j<=n; j++) {	/* use only upper triangle of K,M */
			ri = r[i];
			qj = q[j];
			if ( ri < qj )	{
				Mrq[i][j] = M[ri][qj];
				Drq[i][j] = K[ri][qj] - L*M[ri][qj];
			} else {
				Mrq[i][j] = M[qj][ri];
				Drq[i][j] = K[qj][ri] - L*M[qj][ri];
			}
		}
	}


	/* static condensation of the dynamics matrix */

	xtaiy ( Drq, Drr, Drq, N-n, n, Ac );

	for (i=1; i<=n; i++) {
		for (j=i; j<=n; j++) { /* use only upper triangle of K */
			qi = q[i];
			qj = q[j];
			if ( qi <= qj ) 
			 Kc[j][i]=Kc[i][j] = K[qi][qj]-L*M[qi][qj]-Ac[i][j];
		}
	}


	/* dynamic reduction of the mass matrix */

	xtaiy ( Mrq, Drr, Drq, N-n, n, Ac );

	for (i=1; i<=n; i++) {
		for (j=i; j<=n; j++) { /* use only upper triangle of M */
			qi = q[i];
			qj = q[j];
			if ( qi <= qj )
				Mc[j][i]=Mc[i][j] = M[qi][qj]-Ac[i][j]-Ac[j][i];
		}
	}
	
	aixai ( Drr, Mrr, N-n );

	xAx ( Mrr, Drq, Ac, N-n, n );

	for (i=1; i<=n; i++)	for (j=1; j<=n; j++)	Mc[i][j] += Ac[i][j];

	for (i=1; i<=n; i++)	for (j=1; j<=n; j++)	Kc[i][j] += L*Mc[i][j];

	free_ivector ( r,   1, N-n );
	free_matrix  ( Drr, 1,N-n,1,N-n );
	free_matrix  ( Drq, 1,N-n,1,n );
	free_matrix  ( Mrr, 1,N-n,1,N-n );
	free_matrix  ( Mrq, 1,N-n,1,n );
	free_matrix  ( Ac,  1,n,1,n);

	return;
}


/*---------------------------------------------------------------------------- 
DYN_CONDEN - dynamic condensation of mass and stiffness matrices    8oct01
	     matches the response at a set of frequencies
WARNING: Kc and Mc may be ill-conditioned, and possibly non-positive def.
-----------------------------------------------------------------------------*/
void dyn_conden ( M, K, N, p, n, Mc, Kc, V, f )
float	**M, **K, **Mc, **Kc, **V, *f;
int	N, n, *p;
{
	float	**P, **Pi, **matrix(),
		tmp;
	int	i,j,k, pi, ok; 
	void	pseudo_inv(),	/* pseudo-inverse of a a-symmetric matrix */
		save_matrix(), 
		free_vector(), free_matrix(); 

	P   =  matrix(1,n,1,n);
	Pi  =  matrix(1,n,1,n);

	for (i=1; i<=n; i++)	/* first n modal vectors at primary DoF's */
		for (j=1; j<=n; j++)
			P[i][j] = V[p[i]][j];

	pseudo_inv ( P, Pi, n, n, 1.0e-10 );

	for (i=1; i<=n; i++) { 		/* compute inv(P)' * I * inv(P)	*/
	    for (j=1; j<=n; j++) {
		tmp = 0.0;
	        for (k=1; k<=n; k++)
			tmp += Pi[k][i] * Pi[k][j];
		Mc[i][j] = tmp;
	    }
	}

	for (i=1; i<=n; i++) { 		/* compute inv(P)' * W^2 * inv(P) */
	    for (j=1; j<=n; j++) {
		tmp = 0.0;
	        for (k=1; k<=n; k++) 
			tmp += Pi[k][i] * 4.0*PI*PI*f[k]*f[k] * Pi[k][j];
		Kc[i][j] = tmp;
	    }
	}

	free_matrix  ( P,  1,n,1,n);
	free_matrix  ( Pi, 1,n,1,n);

	return;
}


/*------------------------------------------------------------------------------
XTAIY  -  calculate quadratic form with inverse matrix   X' * inv(A) * Y  
          A is n by n    X is n by m     Y is n by m                    15sep01
------------------------------------------------------------------------------*/
void xtaiy ( X, A, Y, n, m, Ac )
float	**X, **A, **Y, **Ac;
int	n, m;
{
	float	*diag, *b, *x, *vector(), error;
	int	i,j,k, ok, disp=0;
	void	ldl_dcmp(), ldl_mprove(), free_vector();

	diag = vector(1,n);
	x    = vector(1,n);
	b    = vector(1,n);

	for (i=1; i<=n; i++) diag[i] = x[i] = 0.0;

	ldl_dcmp( A, n, diag, b, x, 1, 0, &ok );	/*  L D L'  decomp */

	for (j=1; j<=m; j++) {

		for (k=1; k<=n; k++)  b[k] = Y[k][j];
		ldl_dcmp( A, n, diag, b, x, 0, 1, &ok ); /*  L D L'  bksbtn */

		if (disp)
		 fprintf(stderr,"  RMS matrix error:"); /*improve the solution*/
		error = ok = 1;
		do {
			ldl_mprove ( A, n, diag, b, x, &error, &ok );
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
	free_vector(b,1,n);
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
	float	*diag, *b, *x, **Ai, **XAi, tmp, *vector(), **matrix(), error;
	int	i,j,k, ok, disp=0;
	void	ldl_dcmp(), ldl_mprove(), free_vector(), free_matrix();

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
		 fprintf(stderr,"  RMS matrix error:"); /*improve the solution*/
		error = ok = 1;
		do {
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
			tmp = 0.0;
			for (k=1; k<=n; k++)	tmp += X[i][k]*Ai[k][j];
			XAi[i][j] = tmp;
		}
	}

	for (i=1; i<=n; i++) {		/* compute inv(A) * X * inv(A)	*/
		for (j=1; j<=n; j++) {
			tmp = 0.0;
			for (k=1; k<=n; k++)	tmp += Ai[i][k] * XAi[k][j];
			X[i][j] = tmp;
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
CONTROL_DATA  -  save input data					7nov02
------------------------------------------------------------------------------*/
void control_data ( fp, title, nJ,nM, nF,nD,nR,nW,nP,nT, x,y,z,r, J1,J2, 
		Ax, Asy, Asz, J,Iy,Iz, E,G,p, F,Dp,R,W,P,T, shear, anlyz, geom )
FILE	*fp;
char	*title;
int	nJ, nM, nF, nD, nR, nW, nP, nT, *J1, *J2, *R, shear, anlyz, geom;
float	*x, *y, *z, *r, *Ax,*Asy,*Asz, *J,*Iy,*Iz, *E,*G,*p, *F,*Dp,**W,**P,**T;
{
	int	i,j,n;
        time_t  now;            /* modern time variable type    (DJGPP) */

	(void) time(&now);

	fprintf(fp,"\n");
	for (i=1; i<=80; i++)	fprintf(fp,"_");
	fprintf(fp,"\n-- FRAME version:   1 Mar 2007,");
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
void save_results ( fp, nJ, nM, DoF, J1, J2, F, D, R, Q, err, ok )
FILE	*fp;
int	nJ, nM, DoF, *J1, *J2, *R, ok;
float	*F, *D, **Q, err;
{
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
void modal_results ( fp, nJ, nM, nI, DoF, M, f, V, total_mass, struct_mass,
			iter, sumR, modes, shift, lump, tol, ok)
FILE	*fp;
int	nJ, nM, nI, DoF, iter, sumR, modes, lump, ok;
float	**M, *f, **V, tol, shift, total_mass, struct_mass;
{
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
void mesh ( IO_file, meshfile, plotfile, title, nJ, nM, DoF, x,y,z,L, J1,J2, p, D, exg, anlyz)
char	IO_file[], meshfile[], plotfile[], *title;
int	nJ, nM, DoF, *J1, *J2, anlyz;
float	*x, *y, *z, *L, *p, *D, exg;
{
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
void modal_mesh ( IO_file, meshfile, modefile, plotfile, title, 
			nJ,nM, DoF, modes, x,y,z, L, J1,J2, p, M,f,V, exg,anlyz)
char	IO_file[], meshfile[], modefile[], plotfile[], *title;
int	nJ, nM, DoF, *J1, *J2, modes, anlyz;
float	*x, *y, *z, *L, **M, *f, *p, **V, exg;
{
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
	  fprintf(fpm,"plot '%s' u 2:3 t 'undeformed mesh' w lp ", meshfile );
	  if (!anlyz) fprintf(fpm," lw 2 lt 1 pt 6\n");
	  else fprintf(fpm," lw 1 lt 5 pt 6, '%s' u 1:2 t 'mode-shape %d' w l lw 2 lt 3\n",
								modefl, m );
	  fprintf(fpm,"%c pause -1\n", D3 );
	  fprintf(fpm,"%c set nokey\n", D3 );
	  fprintf(fpm,"%c splot '%s' u 2:3:4 t 'undeformed mesh' w lp ",
								D3, meshfile);
	  if (!anlyz) fprintf(fpm," lw 2 lt 1 pt 6\n");
	  else fprintf(fpm," lw 1 lt 5 pt 6, '%s' u 1:2:3 t 'mode-shape %d' w l lw 2 lt 3\n",
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
void animate( IO_file, meshfile, modefile, plotfile, title, anim,
	nJ,nM, DoF, modes, x,y,z, L, p, J1,J2, f,V, exg, pan )
char	IO_file[], meshfile[], modefile[], plotfile[], *title;
int	nJ, nM, DoF, *J1, *J2, modes, anim[], pan; 
float	*x, *y, *z, *L, *p, *f, **V, exg;
{
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
void bent_beam ( fp, j1, j2, x,y,z, L, p, D, exg )
FILE	*fp;
int	j1, j2;
float	*x, *y, *z, L, p, *D, exg; 
{
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


/*------------------------------------------------------------------------------
DEALLOCATE  -  release the allocated memory				23feb94
------------------------------------------------------------------------------*/
void deallocate( x,y,z,r, L,Le, J1, J2, Ax,Asy,Asz, J,Iy,Iz, E, G,
      K,Q,F,D,R,W,P,T, feF,Fo, d,BMs,JMs,JMx,JMy,JMz, M,f,V, nJ,nM,DoF, modes ) 

int	nJ, nM, DoF, *J1, *J2, *R;
float	*x,*y,*z, *r, *L, *Le, *Ax, *Asy,*Asz, *J,*Iy,*Iz, *E, *G,  **K,
**Q,*F,*D,**W,**P,**T, **feF, *Fo, *d,*BMs,*JMs,*JMx,*JMy,*JMz, **M,*f,**V;
{
	void	free_vector(),	/* memory deallocation for a vector of floats */
		free_matrix(),	/* memory deallocation for a matrix of floats */
		free_ivector();	/* memory deallocation for a vector of ints   */

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


/*------------------------------------------------------------------------------
ITOA  -  Convert an integer n to charcters in s, from K&R, 1978,   p. 59-60
------------------------------------------------------------------------------*/
void	itoa(n,s,k)
int	n,k;
char	s[];
{
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
STRCAT  -  concatenate string t  to string end of string s, K&R, 1978,   p. 44
------------------------------------------------------------------------------*/
int Strcat(s,t)
char	s[], t[];
{
	int	i = 0, j = 0;

	while ( s[i] != '\0' )	i++;	/* find length of s  */
	while ( s[i++] = t[j++] )  ;
}

/*------------------------------------------------------------------------------
STRCPY  -  copy string t  to string s, from K&R, 1978,   p. 101
------------------------------------------------------------------------------*/
int Strcpy(s,t)
char	*s, *t;
{
	while ( *s++ = *t++ )	;
}

/*------------------------------------------------------------------------------
DOTS  -  print a set of dots (periods) 
------------------------------------------------------------------------------*/
void dots(n)
int	n;
{
	int i;
	for (i=1; i<=n; i++)	printf(".");
}
