function [D,R,F,L,Ks] = frame_3dd(XYZ,ELT,RCT,EAIJ,P,U,D)
% [D,R,F,L,Ks] = frame_3dd (XYZ,ELT,RCT,EAIJ,P,U,D)
%
% Solve a a three-dimensional frame analysis problem 
%
% INPUT DATA:
%
%  XYZ : a 4xJ matrix containing the XYZ coordinate of each node
%          row 1 = X-axis coordinate  for each node
%          row 2 = Y-axis coordinate  for each node
%          row 3 = Z-axis coordinate  for each node
%          row 4 = rigid radius       for each node     
%
%  ELT : a 2xB matrix indicating which 2 nodes each frame element connects
%          row 1 = the 'starting' node  for each frame element 
%          row 2 = the 'ending'   node  for each frame element
%
%  RCT : a 6xJ matrix indicated which nodes have reactions
%        0: the node has no reaction in that degree of freedom,
%        1: the node does have a reaction in that degree of freedom.
%
% EAIJ : a 10xB containing the section and modulus properties of each frame el.
%         row 1 = Ax  cross section area                   for each frame el.
%         row 2 = Asy shear area y-direction               for each frame el.
%         row 3 = Asz shear area z-direction               for each frame el. 
%         row 4 = Jxx torsional moment of inertia - x axis for each frame el.
%         row 5 = Iyy bending moment of inertia - y axis   for each frame el.
%         row 6 = Izz bending moment of inertia - z axis   for each frame el.
%         row 7 = E   elastic modulus                      for each frame el.
%         row 8 = G   shear   modulus                      for each frame el.
%         row 9 = p   roll angle                           for each frame el.
%         row 10 = d  mass density                         for each frame el.
%         
%    P : a 6xJ matrix containing the components of the externally applied 
%        forces and moments applied to each node.
%          row 1 = Nodal Force in X-direction    for each node
%          row 2 = Nodal Force in Y-direction    for each node
%          row 3 = Nodal Force in Z-direction    for each node
%          row 4 = Nodal Moment about X-axis     for each node
%          row 5 = Nodal Moment about Y-axis     for each node
%          row 6 = Nodal Moment about Z-axis     for each node
%
%    U : a 3xB matrix containing the unif. dist. load on each frame element
%          row 1 = uniform distributed load along the local element x axis
%          row 2 = uniform distributed load in    the local element y axis
%          row 3 = uniform distributed load in    the local element z axis
%
%    D : a 6xJ matrix of prescribed displacements at the reaction DoF's
%          row 1 = prescribed node displ. in the X-direction for each node
%          row 2 = prescribed node displ. in the Y-direction for each node
%          row 3 = prescribed node displ. in the Z-direction for each node
%          row 4 = prescribed node rot'n  about the X-axis   for each node
%          row 5 = prescribed node rot'n  about the Y-axis   for each node
%          row 6 = prescribed node rot'n  about the Z-axis   for each node
%
% OUTPUT DATA:
%
%    D : a 6xJ matrix   of the deflections and rotations of each node
%    R : a 6xJ matrix   of the reaction forces and moments 
%    F : a 12xB matrix  of the end forces of each frame element 
%    L : a 1xB vector   of the length of each frame element 
%   Ks : a 6Jx6J matrix of the structural stiffness matrix 
%
% http://frame3dd.sourceforge.net/

% NOTE:
% This m-function, frame_3dd.m, executes the system command, frame3dd, 
% to compute the solution.  This m-function interface to frame3dd does not 
% (yet) implement the following features of frame3dd:
%   gravitational loading 
%   point forces applied between the nodes of a member
%   temperature loads
%   multiple load cases
%   modal analysis
%   matrix condensation

% FRAME3DD:
% Static and dynamic structural analysis of 2D and 3D frames and trusses with
% elastic and geometric stiffness.
% ---------------------------------------------------------------------------
% http://www.duke.edu/~hpgavin/frame/
% ---------------------------------------------------------------------------
% Copyright (C) 1992-2008  Henri P. Gavin
% 
%    FRAME3DD is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FRAME3DD is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FRAME3DD.  If not, see <http://www.gnu.org/licenses/>.

% H.P. Gavin, Dept. Civil & Environ. Eng'g, Duke Univ., Mar 31 2009, Apr 8 2011
% *** Please email extensions and enhancements of this function to me. ***

  shear = 0;		     % 1: include shear deformation effects, 0: don't
  geom  = 0;		     % 1: include geometric stiffness effects, 0: don't
  exagg = 10;		     % exaggeration factor for Gnuplot output
  scale = 1.0;		     % zoom scale factor for Gnuplot plotting
  dx = 1.0;		     % x-axis increment for internal force calc's
  IOfilename = 'IOdata';     % name of the Input-Output data file, a text file 
  Mfilename  = 'IOdata_out'; % name of the Matlab results data file, a text file 

  if nargin < 6
         help frame_3dd
         return
  end

  [x,J] = size(XYZ);                    % number of nodes
  [x,B] = size(ELT);                    % number of frame elements

  % --- error checking

  if any(~(size(RCT)==[6,J]))
     error('The dimension of RCT must be 6 by # of nodes.')
  end
  if any(~(size(EAIJ)==[10,B]))
     error('The dimension of EAIJ must be 10 by # of frame elements.')
  end
  if any(~(size(P)==[6,J]))
     error('The dimension of P must be 6 by # of nodes.')
  end
  if any(~(size(U)==[3,B]))
     error('The dimension of U must be 3 by # of frame elements.')
  end
  if any(EAIJ) <= 0
     error('All elements of EAIJ must be greater than zero.')
  end

% open the frame3dd Input-Output Data file
  fp = fopen([IOfilename '.FMM'],'w');

  fprintf(fp,'frame analysis via Matlab interface\n\n');

  fprintf(fp,'%% node data ...\n');
  fprintf(fp,'%4d\t\t%% number of nodes  \n', J);
  fprintf(fp,'%% J\t\tX\t\tY\t\tZ\t\tr\t\tnode coordinates \n');
  for j=1:J
      fprintf(fp,'%4d\t%14.6e\t%14.6e\t%14.6e\t%14.6e\n',j,XYZ(1,j),XYZ(2,j),XYZ(3,j),0);
  end
  fprintf(fp,'\n');

  fprintf(fp,'%% reaction data ...\n');
  nR = sum(max(abs(RCT))~=0);
  fprintf(fp,'%4d    %% number of nodes with reaction forces\n', nR);
  fprintf(fp,'%% j\tRx\tRy\tRz\tRxx\tRyy\tRzz\n'); 
  idx = find(max(abs(RCT)));
  for i=1:nR
      j = idx(i);
      fprintf(fp,'%4d\t%4d\t%4d\t%4d\t%4d\t%4d\t%4d\n', ...
		j, RCT(1,j), RCT(2,j), RCT(3,j), RCT(4,j), RCT(5,j), RCT(6,j) );
  end
  fprintf(fp,'\n');

  fprintf(fp,'%% frame element section property data ...\n');
  fprintf(fp,'%4d\t\t%% number of frame elements\n', B);
  fprintf(fp,'%% m\tn1\tn2\t\tAx\t\tAsy\t\tAsz\t\tJxx\t\tIyy\t\tIzz\t\tE\t\tG\t\tp\tdensity\n');
  for b=1:B
      fprintf(fp,'%4d\t%4d\t%4d\t%14.6e\t%14.6e\t%14.6e\t%14.6e\t%14.6e\t%14.6e\t%14.6e\t%14.6e\t%14.6e\t%14.6e\n', ...
	b, ELT(1,b), ELT(2,b), ...
	EAIJ(1,b), EAIJ(2,b), EAIJ(3,b), ...
	EAIJ(4,b), EAIJ(5,b), EAIJ(6,b), ...
        EAIJ(7,b), EAIJ(8,b), EAIJ(9,b), EAIJ(10,b) );
  end
  fprintf(fp,'\n');

  fprintf(fp,'%4d\t\t%% 1: include shear deformation, 0: do not\n', shear );
  fprintf(fp,'%4d\t\t%% 1: include geometric stiffness, 0: do not\n', geom );
  fprintf(fp,'%14.6e\t%% exagerate deformations in plotting \n', exagg );
  fprintf(fp,'%14.6e\t%% zoom scale factor for 3D plotting \n', scale );
  fprintf(fp,'%14.6e\t%% x-axis increment for internal forces calc\n', dx );
  fprintf(fp,'\n');

  fprintf(fp,'%% static load data ...\n');
  fprintf(fp,'%4d\t\t%% number of static load cases \n', 1);
  fprintf(fp,'\t\t%% begin static load case 1 of 1 \n\n');


  fprintf(fp,'%% gravitational acceleration for self-weight loading\n');
  fprintf(fp,'%% gX         gY         gZ\n');
  fprintf(fp,'  0.0        0.0        0.0\n\n');

  nF = sum(max(abs(P))~=0);
  fprintf(fp,'%4d\t\t%% number of loaded nodes\n', nF);
  fprintf(fp,'%% j\t\tFx\t\tFy\t\tFz\t\tMxx\t\tMyy\t\tMzz\n'); 
  idx = find(max(abs(P)));
  for i=1:nF
      j = idx(i);
      fprintf(fp,'%4d\t%14.6e\t%14.6e\t%14.6e\t%14.6e\t%14.6e\t%14.6e\n', ...
		j, P(1,j), P(2,j), P(3,j), P(4,j), P(5,j), P(6,j) );
  end
  fprintf(fp,'\n');
            
  nU = sum(max(abs(U))~=0);
  fprintf(fp,'%4d\t\t%% number of members with uniform distributed loads \n', nU);
  fprintf(fp,'%% j\t\tUx\t\tUy\t\tUz\n');
  idx = find(max(abs(U)));
  for i=1:nU
      m = idx(i);
      fprintf(fp,'%4d\t%14.6e\t%14.6e\t%14.6e\n', ...
		m, U(1,m), U(2,m), U(3,m) );
  end
  fprintf(fp,'\n');

  nW = 0;
  fprintf(fp,'%4d\t\t%% number of members with trapezoidal loads \n', nW);
  fprintf(fp,'\n');

  nP = 0;
  fprintf(fp,'%4d\t\t%% number of members with internal point loads \n', nP);
  fprintf(fp,'\n');

  nT = 0;
  fprintf(fp,'%4d\t\t%% number of members with temperature loads \n', nT);
  fprintf(fp,'\n');

  nD = sum(max(abs(D))~=0);
  fprintf(fp,'%4d\t\t%% number of nodes with prescribed displacements\n', nD);
  fprintf(fp,'%% j\t\tDx\t\tDy\t\tDz\t\tDxx\t\tDyy\t\tDzz\n'); 
  idx = find(max(abs(D)));
  for i=1:nD
      j = idx(i);
      fprintf(fp,'%4d\t%14.6e\t%14.6e\t%14.6e\t%14.6e\t%14.6e\t%14.6e\n', ...
		j, D(1,j), D(2,j), D(3,j), D(4,j), D(5,j), D(6,j) );
  end
  fprintf(fp,'\n');
  fprintf(fp,'\t\t%% end   static load case 1 of 1 \n\n');

  fprintf(fp,'%% inertial load data ...\n');
  nM = 0;
  fprintf(fp,'%4d\t\t%% number of dynamic modes to analyze \n', nM);
  fprintf(fp,'\n');
 
  fclose(fp);			% close the frame3dd input file

  % compute lengths of each frame element
  L = zeros(1,B);
  for b=1:B
      n1 = ELT(1,b);                          % node 1 of element b
      n2 = ELT(2,b);                          % node 2 of element b

      x1 = XYZ(1,n1); y1 = XYZ(2,n1); z1 = XYZ(3,n1); % coordinates of node 1
      x2 = XYZ(1,n2); y2 = XYZ(2,n2); z2 = XYZ(3,n2); % coordinates of node 2

      L(b) = sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2);
  end

% run the frame3dd analysis
  system(['frame3dd -w -i ' IOfilename '.FMM -o ' IOfilename '.OUT']); 

  run(Mfilename);				% load the results

   D=D1; R=R1; F=F1;				% first load case only

% ------------------------------------------------------------- frame_3dd
