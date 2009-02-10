function [D,R,F,L,Ks] = frame_3dd(XYZ,JTS,RCT,EAIJ,P,W,D)
% [D,R,F,L,Ks] = frame_3dd (XYZ,JTS,RCT,EAIJ,P,W,D)
%
% Solve a a three-dimensional frame analysis problem 
%
% INPUT DATA:
%
%  XYZ : a 4xJ matrix containing the XYZ coordinate of each joint
%          row 1 = X-axis coordinate  for each joint
%          row 2 = Y-axis coordinate  for each joint
%          row 3 = Z-axis coordinate  for each joint
%          row 4 = rigid radius       for each joint     
%
%  JTS : a 2xB matrix indicating which 2 joints each beam element connects
%          row 1 = the 'starting' joint  for each beam
%          row 2 = the 'ending'   joint  for each beam
%
%  RCT : a 6xJ matrix indicated which joints have reactions
%        0: the joint has no reaction in that degree of freedom,
%        1: the joint does have a reaction in that degree of freedom.
%
% EAIJ : a 9xB containing the section and modulus properties of each beam 
%         row 1 = Ax  cross section area                   for each beam
%         row 2 = Asy shear area y-direction               for each beam
%         row 3 = Asz shear area z-direction               for each beam
%         row 4 = Jxx torsional moment of inertia - x axis for each beam
%         row 5 = Iyy bending moment of inertia - y axis   for each beam
%         row 6 = Izz bending moment of inertia - z axis   for each beam
%         row 7 = E   elastic modulus                      for each beam
%         row 8 = G   shear   modulus                      for each beam
%         row 9 = p   roll angle                           for each beam
%         
%    P : a 6xJ matrix containing the components of the external 
%        forces and moments applied to each joint.
%          row 1 = Joint Force in X-direction    for each joint
%          row 2 = Joint Force in Y-direction    for each joint
%          row 3 = Joint Force in Z-direction    for each joint
%          row 4 = Joint Moment about X-axis     for each joint
%          row 5 = Joint Moment about Y-axis     for each joint
%          row 6 = Joint Moment about Z-axis     for each joint
%
%    W : a 3xB matrix containing the unif. dist. load on each beam element
%    D : a 6xJ matrix of prescribed displacements at the reaction DoF's
%
% OUTPUT DATA:
%
%  D : a 6xJ matrix   of the deflections and rotations of each joint
%  R : a 6xJ matrix   of the reaction forces and moments 
%  F : a 12xB matrix  of the end forces of each beam element 
%  L : a 1xB vector   of the length of each beam element 
% Ks : a 6Jx6J matrix of the structural stiffness matrix 
%
% http://www.duke.edu/~hpgavin/frame3dd/

% NOTE:
% This m-function, frame_3dd.m, executes the system command, frame3dd, 
% to compute the solution.  This m-function interface to frame3dd does not 
% (yet) implement the following features of frame3dd:
%   point forces applied between the joints of a member
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

% H.P. Gavin, Dept. Civil & Environ. Eng'g, Duke Univ., September 2008
% *** Please email extensions and enhancements of this function to me. ***

  shear = 0;		     % 1: include shear deformation effects, 0: don't
  geom  = 0;		     % 1: include geometric stiffness effects, 0: don't
  exagg = 10;		     % exaggeration factor for Gnuplot output
  anlyz = 1;		     % 1: run a frame3dd analysis, 0: data check only
  IOfilename = 'IOdata';     % name of the Input-Output data file, a text file 
  Mfilename  = 'IOdata_out'; % name of the Matlab results data file, a text file 

  if nargin < 6
         help frame_3dd
         return
  end

  [x,J] = size(XYZ);                    % number of joints
  [x,B] = size(JTS);                    % number of beam elements

  % --- error checking

  if any(~(size(RCT)==[6,J]))
     error('The dimension of RCT must be 6 by # of joints.')
  end
  if any(~(size(EAIJ)==[9,B]))
     error('The dimension of EAIJ must be 9 by # of beams.')
  end
  if any(~(size(P)==[6,J]))
     error('The dimension of P must be 6 by # of joints.')
  end
  if any(~(size(W)==[3,B]))
     error('The dimension of W must be 3 by # of beams.')
  end
  if any(EAIJ) <= 0
     error('All elements of EAIJ must be greater than zero.')
  end

% open the frame3dd Input-Output Data file
  fp = fopen([IOfilename '.FMM'],'w');

  fprintf(fp,'frame analysis\n\n');

  fprintf(fp,'%% joint data ...\n');
  fprintf(fp,'%d\t\t%% number of joints  \n', J);
  fprintf(fp,'%% J\t\tX\t\tY\t\tZ\t\tr\t\tjoint coordinates \n');
  for j=1:J
      fprintf(fp,'%d\t%e\t%e\t%e\t%e\n',j,XYZ(1,j),XYZ(2,j),XYZ(3,j),0);
  end
  fprintf(fp,'\n');

  fprintf(fp,'%% reaction data ...\n');
  nR = sum(max(abs(RCT))~=0);
  fprintf(fp,'%d    %% number of joints with reaction forces\n', nR);
  fprintf(fp,'%% j\tRx\tRy\tRz\tRxx\tRyy\tRzz\n'); 
  idx = find(max(abs(RCT)));
  for i=1:nR
      j = idx(i);
      fprintf(fp,'%d\t%d\t%d\t%d\t%d\t%d\t%d\n', ...
		j, RCT(1,j), RCT(2,j), RCT(3,j), RCT(4,j), RCT(5,j), RCT(6,j) );
  end
  fprintf(fp,'\n');

  fprintf(fp,'%% beam section property data ...\n');
  fprintf(fp,'%d\t\t%% number of beam elements\n', B);
  fprintf(fp,'%% m\tj1\tj2\t\tAx\t\tAsy\t\tAsz\t\tJxx\t\tIyy\t\tIzz\t\tE\t\tG\t\tp\n');
  for b=1:B
      fprintf(fp,'%d\t%d\t%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n', ...
	b, JTS(1,b), JTS(2,b), ...
	EAIJ(1,b), EAIJ(2,b), EAIJ(3,b), ...
	EAIJ(4,b), EAIJ(5,b), EAIJ(6,b), ...
        EAIJ(7,b), EAIJ(8,b), EAIJ(9,b) );
  end
  fprintf(fp,'\n');

  fprintf(fp,'%d\t\t%% 1: include shear deformation, 0: do not\n', shear );
  fprintf(fp,'%d\t\t%% 1: include geometric stiffness, 0: do not\n', geom );
  fprintf(fp,'%e\t%% exxagerate deformations in plotting \n', exagg );
  fprintf(fp,'%d\t\t%% 1: frame3dd analysis, 0: data check only\n', anlyz);
  fprintf(fp,'\n');

  fprintf(fp,'%% static load data ...\n');
  fprintf(fp,'%d\t\t%% number of static load cases \n', 1);
  fprintf(fp,'\t\t%% begin static load case 1 of 1 \n');
  nF = sum(max(abs(P))~=0);
  fprintf(fp,'%d\t\t%% number of loaded joints\n', nF);
  fprintf(fp,'%% j\t\tFx\t\tFy\t\tFz\t\tMxx\t\tMyy\t\tMzz\n'); 
  idx = find(max(abs(P)));
  for i=1:nF
      j = idx(i);
      fprintf(fp,'%d\t%e\t%e\t%e\t%e\t%e\t%e\n', ...
		j, P(1,j), P(2,j), P(3,j), P(4,j), P(5,j), P(6,j) );
  end
  fprintf(fp,'\n');
            
  nW = sum(max(abs(W))~=0);
  fprintf(fp,'%d\t\t%% number of members with distributed loads \n', nW);
  fprintf(fp,'%% j\t\tWx\t\tWy\t\tWz\n');
  idx = find(max(abs(W)));
  for i=1:nW
      m = idx(i);
      fprintf(fp,'%d\t%e\t%e\t%e\n', ...
		m, W(1,m), W(2,m), W(3,m) );
  end
  fprintf(fp,'\n');

  nP = 0;
  fprintf(fp,'%d\t\t%% number of members with internal point loads \n', nP);
  fprintf(fp,'\n');

  nT = 0;
  fprintf(fp,'%d\t\t%% number of members with temperature loads \n', nT);
  fprintf(fp,'\n');

  nD = sum(max(abs(D))~=0);
  fprintf(fp,'%d\t\t%% number of joints with prescribed displacements\n', nD);
  fprintf(fp,'%% j\t\tDx\t\tDy\t\tDz\t\tDxx\t\tDyy\t\tDzz\n'); 
  idx = find(max(abs(D)));
  for i=1:nD
      j = idx(i);
      fprintf(fp,'%d\t%e\t%e\t%e\t%e\t%e\t%e\n', ...
		j, D(1,j), D(2,j), D(3,j), D(4,j), D(5,j), D(6,j) );
  end
  fprintf(fp,'\n');
  fprintf(fp,'\t\t%% end   static load case 1 of 1 \n\n');

  fprintf(fp,'%% inertial load data ...\n');
  nM = 0;
  fprintf(fp,'%d\t\t%% number of dynamic modes to analyze \n', nM);
  fprintf(fp,'\n');
 
  fclose(fp);			% close the frame3dd input file

  % compute lengths of each beam element
  L = zeros(1,B);
  for b=1:B
      j1 = JTS(1,b);                          % joint 1 of bar b
      j2 = JTS(2,b);                          % joint 2 of bar b

      x1 = XYZ(1,j1); y1 = XYZ(2,j1); z1 = XYZ(3,j1); % coordinates of joint 1
      x2 = XYZ(1,j2); y2 = XYZ(2,j2); z2 = XYZ(3,j2); % coordinates of joint 2

      L(b) = sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2);
  end

  system(['frame3dd ' IOfilename '.FMM']); 	% run the frame3dd analysis

  run(Mfilename);				% load the results

   D=D1; R=R1; F=F1;				% first load case only

% ------------------------------------------------------------- frame_3dd
