Template for input file format for FRAME program

nJ						% number of joints
nM						% number of members
nL						% number of load cases

  J[1]	x[1]	y[1]	z[1]	r[1]	
  :	:	:	:	:		 %joint numbers and coordinates
  J[nJ]	x[nJ]	y[nJ]	z[nJ]	r[1]
					    %member numbers, loc'ns, and prop's
  M[1]   J1[1]  J2[1]   Ax[1]  Asy[1]  Asz[1]  Jp[1]  Iy[1]  Iz[1]   E[1]  G[1]  p[1]
  :      :      :       :       :       :       :      :      :       :     :     :
  M[nM]  J1[nM] J2[nM]  Ax[nM] Asy[nM] Asz[nM] Jp[nM] Iy[nM] Iz[nM]  E[nM] G[nM] p[nM]

nR					       %number of joints with reactions
  J[1]	Rx[1]	Ry[1]	Rz[1]	Rxx[1]	Ryy[1]	Rzz[1]
  :	:	:	:	:	:	:	       %0:free, 1:fixed
  J[nR]	Rx[nR]	Ry[nR]	Rz[nR]	Rxx[nR]	Ryy[nR]	Rzz[nR]


shear						  %1: include shear deformation
mesh_file						   %mesh data file name
ann_file					     %mesh annotation file name
exagg						  %exaggerate mesh deformations
anlyz				     %1: stiffness analysis, 0: data check only


                                                       % begin load case 1 of nL
nF						       %number of loaded joints
  J[1]	Fx[1]	Fy[1]	Fz[1]	Mxx[1]	Myy[1]	Mzz[1]
  :	:	:	:	:	:	:		  %nodal forces
  J[nF]	Fx[nF]	Fy[nF]	Fz[nF]	Mxx[nF]	Myy[nF]	Mzz[nF]

nW					   %number of uniform distributed loads
  M[1]	Wx[1]	Wy[1]	Wz[1]
  :	:	:	:	    %uniform member loads in member coordinates
  M[nW]	Wx[nW]	Wy[nW]	Wz[nW]

nP					    %number of concentrated point loads
  M[1]	Px[1]	Py[1]	Pz[1]	x[1]	    %point loads in member coordinates 
  :	:	:	:	:	    %and x=distance from coordinate J1 
  M[nP]	Px[nP]	Py[nP]	Pz[nP]	x[nP]

nT					 %number of members temperature changes
  M[1]	a[1]	hy[1]	hz[1]	Ty+[1]	Ty-[1]	Tz+[1]	Tz-[1]	%member no.,   
  :	:	:	:	:	:	:	:	%temp. coef.
  M[nT]	a[nT]	hy[nT]	hz[nT]	Ty+[nT]	Ty-[nT]	Tz+[nT]	Tz-[nT] %sizes, & temps

nD		         %number of joints with prescribed displacements nD<=nR
  J[1]	Dx[1]	Dy[1]	Dz[3]	Dxx[1]	Dyy[1]	Dzz[1]
  :	:	:	:	:	:	:     %prescribed displacements
  J[nD]	Dx[nD]	Dy[nD]	Dz[nD]	Dxx[nD]	Dyy[nD]	Dzz[nD]
                                                       % end load case 1 of nL

                                                       % begin load case 2 of nL
nF						       %number of loaded joints
  J[1]	Fx[1]	Fy[1]	Fz[1]	Mxx[1]	Myy[1]	Mzz[1]
  :	:	:	:	:	:	:		  %nodal forces
  J[nF]	Fx[nF]	Fy[nF]	Fz[nF]	Mxx[nF]	Myy[nF]	Mzz[nF]

nW					   %number of uniform distributed loads
  M[1]	Wx[1]	Wy[1]	Wz[1]
  :	:	:	:	    %uniform member loads in member coordinates
  M[nW]	Wx[nW]	Wy[nW]	Wz[nW]

nP					    %number of concentrated point loads
  M[1]	Px[1]	Py[1]	Pz[1]	x[1]	    %point loads in member coordinates 
  :	:	:	:	:	    %and x=distance from coordinate J1 
  M[nP]	Px[nP]	Py[nP]	Pz[nP]	x[nP]

nT					 %number of members temperature changes
  M[1]	a[1]	hy[1]	hz[1]	Ty+[1]	Ty-[1]	Tz+[1]	Tz-[1]	%member no.,   
  :	:	:	:	:	:	:	:	%temp. coef.
  M[nT]	a[nT]	hy[nT]	hz[nT]	Ty+[nT]	Ty-[nT]	Tz+[nT]	Tz-[nT] %sizes, & temps

nD		         %number of joints with prescribed displacements nD<=nR
  J[1]	Dx[1]	Dy[1]	Dz[3]	Dxx[1]	Dyy[1]	Dzz[1]
  :	:	:	:	:	:	:     %prescribed displacements
  J[nD]	Dx[nD]	Dy[nD]	Dz[nD]	Dxx[nD]	Dyy[nD]	Dzz[nD]
                                                       % end load case 2 of nL

% repeat load case information for each load case


modes						       %number of desired modes
lump						 %0: consistent mass, 1: lumped
mode_file					     %mode shape data file name
tol						  %convergence tolerance ~ 1e-4
shift		i   %shift-factor for rigid body modes, make 0 for pos.def. [K]
  M[1]	d[1]	BMs[1]
  :	:	:  %beam density and extra beam masses, not including self mass
  M[nM]	d[nM]	BMs[nM]

  J[1]	JMs[1]    JMx[1]   JMy[1]   JMz[1]  %joint masses and rotatory inertias
  :	:	  :	   :	    :			    %global coordinates
  J[nJ]	JMs[nJM]  JMx[nJ]  JMy[nJ]  JMz[nJ]

anim[0] ... anim[m]   %list of modes to be animated, sorted by increasing freq.


