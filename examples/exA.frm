Example A: linear static analysis of a 2D truss with support settlement (kips,in) 

12				% number of joints 
%joint  x       y       z       r                 units: inches

 1	0.0	0.0	0.0	0.0
 2	12.0	0.0	0.0	0.0
 3	24.0	0.0	0.0	0.0
 4	36.0	0.0	0.0	0.0
 5	48.0	0.0	0.0	0.0
 6	60.0	0.0	0.0	0.0
 7	72.0	0.0	0.0	0.0
 8	12.0	12.0	0.0	0.0
 9	24.0	12.0	0.0	0.0
10	36.0	12.0	0.0	0.0
11	48.0	12.0	0.0	0.0
12	60.0	12.0	0.0	0.0

12				% number of joints with reactions
%J     x y z xx yy zz		1=fixed, 0=free

  1	1 1 1  1  1  0
  2	0 0 1  1  1  0
  3	0 0 1  1  1  0 
  4	0 0 1  1  1  0
  5	0 0 1  1  1  0
  6	0 0 1  1  1  0
  7	0 1 1  1  1  0
  8	1 0 1  1  1  0
  9	0 0 1  1  1  0
 10	0 0 1  1  1  0
 11	0 0 1  1  1  0
 12	0 0 1  1  1  0

21				% number of members 
%m j1 j2 Ax    Asy     Asz     Jxx     Iyy     Izz       E      G     p
%        in^2  in^2    in^2    in^4    in^4    in^4      ksi    ksi   deg. 

 1 1 2	30.0	1.0	1.0	1.0	1.0	0.01	29000 11500    0
 2 2 3	30.0	1.0	1.0	1.0	1.0	0.01	29000 11500    0
 3 3 4	30.0	1.0	1.0	1.0	1.0	0.01	29000 11500    0
 4 4 5	30.0	1.0	1.0	1.0	1.0	0.01	29000 11500    0
 5 5 6	30.0	1.0	1.0	1.0	1.0	0.01	29000 11500    0
 6 6 7	30.0	1.0	1.0	1.0	1.0	0.01	29000 11500    0
 7 1 8	30.0	1.0	1.0	1.0	1.0	0.01	29000 11500    0
 8 2 8	30.0	1.0	1.0	1.0	1.0	0.01	29000 11500    0
 9 2 9	30.0	1.0	1.0	1.0	1.0	0.01	29000 11500    0
10 3 9	30.0	1.0	1.0	1.0	1.0	0.01	29000 11500    0
11 4 9	30.0	1.0	1.0	1.0	1.0	0.01	29000 11500    0
12 4 10	30.0	1.0	1.0	1.0	1.0	0.01	29000 11500    0
13 4 11	30.0	1.0	1.0	1.0	1.0	0.01	29000 11500    0
14 5 11	30.0	1.0	1.0	1.0	1.0	0.01	29000 11500    0
15 6 11	30.0	1.0	1.0	1.0	1.0	0.01	29000 11500    0
16 6 12	30.0	1.0	1.0	1.0	1.0	0.01	29000 11500    0
17 7 12	30.0	1.0	1.0	1.0	1.0	0.01	29000 11500    0
18 8  9	30.0	1.0	1.0	1.0	1.0	0.01	29000 11500    0
19 9 10	30.0	1.0	1.0	1.0	1.0	0.01	29000 11500    0
20 10 11 30.0	1.0	1.0	1.0	1.0	0.01	29000 11500    0
21 11 12 30.0	1.0	1.0	1.0	1.0	0.01	29000 11500    0

 
0                               % 1: include shear deformation
0                               % 1: include geometric stiffness
10.0                            % exaggerate mesh deformations
1                               % 1: stiffness analysis, 0: data check only


2				% number of static load cases
				% Begin Static Load Case 1 of 2
5				% number of loaded joints
%J      Fx       Fy     Fz      Mxx     Myy     Mzz
%       kip      kip    kip     in.k    in.k    in.k
 2	0.0	-20.0	0.0	0.0	0.0	0.0
 3	0.0	-40.0	0.0	0.0	0.0	0.0
 4	0.0	-30.0  	0.0	0.0	0.0	0.0
 5	0.0	-10.0  	0.0	0.0	0.0	0.0
 6	0.0	-30.0  	0.0	0.0	0.0	0.0

0				% number of uniform loads
0				% number of trapezoidal loads
0				% number of internal concentrated loads
0				% number of temperature loads

1				% number of support settlements
%J     Dx      Dy      Dz      Dxx     Dyy     Dzz
%      in      in      in      rad.    rad.    rad.
  8 	0.01	0.0	0.0	0.0	0.0	0.0
				% End   Static Load Case 1 of 2

				% Begin Static Load Case 2 of 2
3				% number of loaded joints
%J      Fx       Fy     Fz      Mxx     Myy     Mzz
%       kip      kip    kip     in.k    in.k    in.k
 3	20.0	 0.0	0.0	0.0	0.0	0.0
 4	30.0	 0.0  	0.0	0.0	0.0	0.0
 5	20.0	 0.0  	0.0	0.0	0.0	0.0

0				% number of uniform loads
0				% number of trapezoidal loads
0				% number of internal concentrated loads

3				% number of temperature loads
%M  a      hy   hz   Ty+  Ty-  Tz+  Tz- 
%   /.F    in   in   F    F    F    F
10 6e-12   1.5  1.5  10	  10   10   10
13 6e-12   1.5  1.5  15   15   15   15  
15 6e-12   1.5  1.5  17   17   17   17  

2				% number of support settlements
%J     Dx      Dy      Dz      Dxx     Dyy     Dzz
%      in      in      in      rad     rad     rad
  1 	0.0    -0.01	0.0	0.0	0.0	0.0
  8 	0.05	0.0	0.0	0.0	0.0	0.0
				% End   Static Load Case 2 of 2


0				% number of desired dynamic modes of vibration


