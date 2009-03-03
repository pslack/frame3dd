Example I: a triangular tower  (kip,in)

15			% number of joints 

%.joint    x         y         z         r
%          in        in        in        in
 
 1      -100         0         0         0
 2       100         0         0         0
 3         0        70         0         0         
 4      -100         0        80         0
 5       100         0        80         0
 6         0        70        80         0         
 7      -100         0       180         0
 8       100         0       180         0
 9         0        70       180         0         
10      -100         0       310         0
11       100         0       310         0
12         0        70       310         0         
13      -100         0       510         0
14       100         0       510         0
15         0        70       510         0         

3                               % number of joints with reactions
%.J     x  y  z xx yy zz          1=fixed, 0=free
  1     1  1  1  1  1  1
  2     1  1  1  1  1  1
  3     1  1  1  1  1  1


24			% number of members
%.m    j1 j2   Ax  Asy   Asz   Jxx    Iyy   Izz    E      G    p
%             in^2 in^2  in^2  in^2   in^2  in^2   ksi    ksi  deg

 1      1  4  100  60    60    500   1000  1000    9990   3800 0
 2      2  5  100  60    60    500   1000  1000    9990   3800 0
 3      3  6  100  60    60    500   1000  1000    9990   3800 0
 4      4  7  100  60    60    500   1000  1000    9990   3800 0
 5      5  8  100  60    60    500   1000  1000    9990   3800 0
 6      6  9  100  60    60    500   1000  1000    9990   3800 0
 7      7 10  100  60    60    500   1000  1000    9990   3800 0
 8      8 11  100  60    60    500   1000  1000    9990   3800 0
 9      9 12  100  60    60    500   1000  1000    9990   3800 0
10     10 13  100  60    60    500   1000  1000    9990   3800 0
11     11 14  100  60    60    500   1000  1000    9990   3800 0
12     12 15  100  60    60    500   1000  1000    9990   3800 0

13      4  5  100  60    60    500   1000  1000    9990   3800 0
14      5  6  100  60    60    500   1000  1000    9990   3800 0
15      6  4  100  60    60    500   1000  1000    9990   3800 0
16      7  8  100  60    60    500   1000  1000    9990   3800 0
17      8  9  100  60    60    500   1000  1000    9990   3800 0
18      9  7  100  60    60    500   1000  1000    9990   3800 0
19     10 11  100  60    60    500   1000  1000    9990   3800 0
20     11 12  100  60    60    500   1000  1000    9990   3800 0
21     12 10  100  60    60    500   1000  1000    9990   3800 0
22     13 14  100  60    60    500   1000  1000    9990   3800 0
23     14 15  100  60    60    500   1000  1000    9990   3800 0
24     15 13  100  60    60    500   1000  1000    9990   3800 0

1                               % 1: include shear deformation
1                               % 1: include geometric stiffness
2.0                             % exaggerate mesh deformations
1                               % 1: stiffness analysis, 0: data check only


1			% number of static load cases
				% Begin Static Load Case 1 of 1
1                               % number of loaded joints
%.J       Fx         Fy     Fz     Mxx     Myy     Mzz
%          k         k      k       k.in   k.in    k.in
  15       0        200      0       0       0       0

12                               % number of uniform distributed loads
%..j      Ux         Uy      Uz
%                            k/in
  13      0          0      -2.361
  14      0          0      -2.361
  15      0          0      -2.361
  16      0          0      -2.361
  17      0          0      -2.361
  18      0          0      -2.361
  19      0          0      -2.361
  20      0          0      -2.361
  21      0          0      -2.361
  22      0          0      -2.361
  23      0          0      -2.361
  24      0          0      -2.361

0                               % number of trapezoidal distributed loads
0                               % number of internal concentrated loads
0                               % number of temperature loads
0                               % number of support settlements
				% End   Static Load Case 1 of 1


4                               % number of desired dynamic modes of vibration
1                               % 1: subspace Jacobi     2: Stodola
0                               % 0: consistent mass ... 1: lumped mass matrix
1e-6                            % mode shape tolerance
0.0                             % shift value ... for unrestrained structures

%.M    density mass
%      k/in^3  kip
 1     2.52e-7  0			% bar numbers, density, and extra mass
 2     2.52e-7  0
 3     2.52e-7  0
 4     2.52e-7  0
 5     2.52e-7  0
 6     2.52e-7  0
 7     2.52e-7  0
 8     2.52e-7  0
 9     2.52e-7  0
10     2.52e-7  0
11     2.52e-7  0
12     2.52e-7  0
13     2.52e-7  0
14     2.52e-7  0
15     2.52e-7  0
16     2.52e-7  0
17     2.52e-7  0
18     2.52e-7  0
19     2.52e-7  0
20     2.52e-7  0
21     2.52e-7  0
22     2.52e-7  0
23     2.52e-7  0
24     2.52e-7  0

0                % number of joints with extra mass or inertia
%.j    M        Ixx   Iyy    Izz -  joints and concentrated mass and inertia

4                               % number of modes to animate
 1  2  3  4                     % modes to animate
1                               % pan rate during animation


