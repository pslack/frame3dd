Example E: a three dimensional structure showing lateral-torsional dynamic modes

12			% number of joints
%.joint  x       y       z       r


   1     0       0       0       0
   2    72       0       0       0
   3   144       0       0       0
   4   144      36       0       0
   5   144      72       0       0
   6    72      72       0       0
   7     0      72       0       0
   8     0      36       0       0
   9     0       0    -120       0
  10   144       0    -120       0
  11    72      72    -120       0
  12    72      36       0       0

3                               % number of joints with reactions
%.J     x    y    z   xx   yy   zz          1=fixed, 0=free
  9     1    1    1    1    1    1
 10     1    1    1    1    1    1
 11     1    1    1    1    1    1

15			% number of members
%.m     j1     j2      Ax      Asy   Asz      Jxx    Iyy    Izz    E      G     p

  1      1      2      1100    800   800      1000    500    500   999000 11100 0
  2      2      3      1100    800   800      1000    500    500   999000 11100 0
  3      3      4      1100    800   800      1000    500    500   999000 11100 0
  4      4      5      1100    800   800      1000    500    500   999000 11100 0
  5      5      6      1100    800   800      1000    500    500   999000 11100 0
  6      6      7      1100    800   800      1000    500    500   999000 11100 0
  7      7      8      1100    800   800      1000    500    500   999000 11100 0
  8      8      1      1100    800   800      1000    500    500   999000 11100 0
  9      9      1      1100    800   800       001    110    110    29000 11100 0
 10     10      3      1100    800   800       001    110    110    29000 11100 0
 11     11      6      1100    800   800       001    110    110    29000 11100 0
 12     12      2      1100    800   800      1000    500    500   999000 11100 0
 13     12      4      1100    800   800      1000    500    500   999000 11100 0
 14     12      6      1100    800   800      1000    500    500   999000 11100 0
 15     12      8      1100    800   800      1000    500    500   999000 11100 0
  

1                               % 1: include shear deformation
1                               % 1: include geometric stiffness
20.0                            % exaggerate mesh deformations
1                               % 1: stiffness analysis, 0: data check only


1			% number of static load cases
				% Begin Static Load Case 1 of 1
1                               % number of loaded joints
%..J       Fx Fy   Fz   Mxx  Myy  Mzz
   3        0  1  -10   0    0    0

0                               % number of distributed loads
0                               % number of internal concentrated loads
0                               % number of members with temperature loads
0                               % number of joints with support settlements
				% End   Static Load Case 1 of 1


4                               % number of desired dynamic modes of vibration
1                               % 1: subspace Jacobi     2: Stodola
0                               % 0: consistent mass ... 1: lumped mass matrix
1e-5                            % mode shape tolerance
1.0                             % shift value ... for unrestrained structures


  1     7.32e-19    0		% bar numbers, density, and extra mass
  2     7.32e-19    0
  3     7.32e-19    0
  4     7.32e-19    0
  5     7.32e-19    0
  6     7.32e-19    0
  7     7.32e-19    0
  8     7.32e-19    0
  9     7.32e-19    0
 10     7.32e-19    0
 11     7.32e-19    0
 12     7.32e-19    0
 13     7.32e-19    0
 14     7.32e-19    0
 15     7.32e-19    0

1                                 % number of joints with extra mass or inertia
%.j     M       Ixx      Iyy      Izz - joints and concentrated mass and inertia
 12     3.388     0        0        839.37


4                               % number of modes to animate
 1  2  3  4                     % modes to animate
0                               % don't pan during animation

2    % Condensation Method:   0= none   1= static   2= Guyan   3= Dynamic
1                               % number of condensed joints
  12    1  1  0   0  0  1	% joint number, 1: condense dof, 0: don't 

  1 2 3 			% modes to match for dynamic condensation

