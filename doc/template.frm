Template for input file format for FRAME program

% joint data ...
nJ						% number of joints
  J[1]	x[1]	y[1]	z[1]	r[1]	
  :	:	:	:	:		 % joint numbers and coordinates
  J[nJ]	x[nJ]	y[nJ]	z[nJ]	r[1]

% beam element property data ...
nB						% number of beam elements 
					    % member numbers, loc'ns, and prop's
  M[1]   J1[1]  J2[1]   Ax[1]  Asy[1]  Asz[1]  Jp[1]  Iy[1]  Iz[1]   E[1]  G[1]  p[1]
  :      :      :       :       :       :       :      :      :       :     :     :
  M[nM]  J1[nM] J2[nM]  Ax[nM] Asy[nM] Asz[nM] Jp[nM] Iy[nM] Iz[nM]  E[nM] G[nM] p[nM]

% reaction data ...
nR					       % number of joints with reactions
  J[1]	Rx[1]	Ry[1]	Rz[1]	Rxx[1]	Ryy[1]	Rzz[1]
  :	:	:	:	:	:	:	       % 0:free, 1:fixed
  J[nR]	Rx[nR]	Ry[nR]	Rz[nR]	Rxx[nR]	Ryy[nR]	Rzz[nR]

shear						  % 1: include shear deformation
/tmp/mesh_file						   % mesh data file name
plot_file					     % mesh annotation file name
exagg						  % exaggerate mesh deformations
anlyz				     % 1: stiffness analysis, 0: data check only

% load data ...
nL						% number of static load cases
% Begin Static Load Case 1
nF						       % number of loaded joints
  J[1]	Fx[1]	Fy[1]	Fz[1]	Mxx[1]	Myy[1]	Mzz[1]
  :	:	:	:	:	:	:		  % nodal forces
  J[nF]	Fx[nF]	Fy[nF]	Fz[nF]	Mxx[nF]	Myy[nF]	Mzz[nF]

nW					   % number of uniform distributed loads
  M[1]	Wx[1]	Wy[1]	Wz[1]
  :	:	:	:	    % uniform member loads in member coordinates
  M[nW]	Wx[nW]	Wy[nW]	Wz[nW]

nP					    % number of concentrated point loads
  M[1]	Px[1]	Py[1]	Pz[1]	x[1]	    % point loads in member coordinates 
  :	:	:	:	:	    % and x=distance from coordinate J1 
  M[nP]	Px[nP]	Py[nP]	Pz[nP]	x[nP]

nT					 % number of members temperature changes
  M[1]	a[1]	hy[1]	hz[1]	Ty+[1]	Ty-[1]	Tz+[1]	Tz-[1]	% member no.,   
  :	:	:	:	:	:	:	:	% temp. coef.
  M[nT]	a[nT]	hy[nT]	hz[nT]	Ty+[nT]	Ty-[nT]	Tz+[nT]	Tz-[nT] % sizes, & temps

nD		         % number of joints with prescribed displacements nD<=nR
  J[1]	Dx[1]	Dy[1]	Dz[3]	Dxx[1]	Dyy[1]	Dzz[1]
  :	:	:	:	:	:	:     % prescribed displacements
  J[nD]	Dx[nD]	Dy[nD]	Dz[nD]	Dxx[nD]	Dyy[nD]	Dzz[nD]
% End Static Load Case 1

% Begin Static Load Case 2
nF						       % number of loaded joints
  J[1]	Fx[1]	Fy[1]	Fz[1]	Mxx[1]	Myy[1]	Mzz[1]
  :	:	:	:	:	:	:		  % nodal forces
  J[nF]	Fx[nF]	Fy[nF]	Fz[nF]	Mxx[nF]	Myy[nF]	Mzz[nF]

nW					   % number of uniform distributed loads
  M[1]	Wx[1]	Wy[1]	Wz[1]
  :	:	:	:	    % uniform member loads in member coordinates
  M[nW]	Wx[nW]	Wy[nW]	Wz[nW]

nP					    % number of concentrated point loads
  M[1]	Px[1]	Py[1]	Pz[1]	x[1]	    % point loads in member coordinates 
  :	:	:	:	:	    % and x=distance from coordinate J1 
  M[nP]	Px[nP]	Py[nP]	Pz[nP]	x[nP]

nT					 % number of members temperature changes
  M[1]	a[1]	hy[1]	hz[1]	Ty+[1]	Ty-[1]	Tz+[1]	Tz-[1]	% member no.,   
  :	:	:	:	:	:	:	:	% temp. coef.
  M[nT]	a[nT]	hy[nT]	hz[nT]	Ty+[nT]	Ty-[nT]	Tz+[nT]	Tz-[nT] % sizes, & temps

nD		         % number of joints with prescribed displacements nD<=nR
  J[1]	Dx[1]	Dy[1]	Dz[3]	Dxx[1]	Dyy[1]	Dzz[1]
  :	:	:	:	:	:	:     % prescribed displacements
  J[nD]	Dx[nD]	Dy[nD]	Dz[nD]	Dxx[nD]	Dyy[nD]	Dzz[nD]
% End Static Load Case 2

% repeat up to 128 static load cases

% dynamic analysis data ...
modes					              % number of desired modes
Mmethod                                        % 1: Subspace Jacobi, 2: Stodola
lump					        % 0: consistent mass, 1: lumped
/tmp/mode_file					    % mode shape data file name
tol		  			         % convergence tolerance ~ 1e-4
shift		   % shift-factor for rigid body modes, make 0 for pos.def. [K]

 M[1]  d[1]  BMs[1]
   :     :       :                           % beam density and extra beam mass
 M[nM] d[nM] BMs[nM]

nI                          % number of joints with extra joint mass or inertia
 J[1]  JMs[1]   JMx[1]  JMy[1]  JMz[1]     % joint masses and rotatory inertias 
   :       :        :       :       :                      % global coordinates 
 J[nI] JMs[nI]  JMx[nI] JMy[nI] JMz[nI]

% animation data ...
nA						% number of modes to be animated
 anim[0] ... anim[nA] % list of modes to be animated, sorted by increasing freq 
pan                                         % 1: pan during animation; 0: don't

% matrix condensation data ...
Cmethod % matrix condensation method ...  0=none, 1=static, 2=Guyan, 3=dynamic
nC                                                % number of condensed joints
 J[1]  cx[1]  cy[1]  cz[1]   cxx[1]  cyy[1] czz[1]
    :      :      :      :        :       :       :    % 1: condense; 0: don't
 J[nC] cx[nC] cy[nC] cz[nC]  cxx[nC] cyy[nC]  czz[nC]

  m[1]   m[2]   m[3]  ...      % list of modes matched in dynamic condensation
                               % if Cmethod == 1, only mode m[1] is matched.

