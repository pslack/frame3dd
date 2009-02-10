Example B: a pyramid-shaped frame --- static and dynamic analysis

5				% number of joints 
%.joint  x       y       z       r

1	0.0	0.0	100.0	0.0
2	-120.0	-90.0	0.0	0.0	
3	 120.0	-90.0	0.0	0.0	
4	 120.0	 90.0	0.0	0.0	
5	-120.0	 90.0	0.0	0.0	

4                               % number of joints with reactions
%.J     x y z xx yy zz          1=fixed, 0=free

  2	1 1 1 1 1 1
  3	1 1 1 1 1 1
  4	1 1 1 1 1 1
  5	1 1 1 1 1 1

4				% number of members			
%.m j1 j2 Ax    Asy     Asz     Jxx     Iyy     Izz         E      G   p

1 1 2	10.	8.0	8.0	500.0	300.0	200.0	1000.0	700.0  0
2 1 3	10.	8.0	8.0	500.0	300.0	200.0	1000.0	700.0  0
3 1 4	10.	8.0	8.0	500.0	300.0	200.0	1000.0	700.0  0
4 1 5	10.	8.0	8.0	500.0	300.0	200.0	1000.0	700.0  0

 
1                               % 1: include shear deformation
1                               % 1: include geometric stiffness
10.0                            % exaggerate mesh deformations
1                               % 1: stiffness analysis, 0: data check only

3				% number of static load cases
				% Begin Static Load Case 1 of 3
1				% number of loaded joints
%.J      Fx       Fy     Fz      Mxx     Myy     Mzz
 1	10.00	-20.00	-100	0.0	0.0	0.0
0                               % number of distributed loads
0                               % number of internal concentrated loads
0                               % number of members with temperature loads
0                               % number of joints with support settlements
				% End   Static Load Case 1 of 3

				% Begin Static Load Case 2 of 3
0				% number of loaded joints
2                               % number of distributed loads
%.M    Wx   Wy   Wz
  2    0    0.1  0
  1    0    0    0.1
0                               % number of internal concentrated loads
1                               % number of members with temperature loads
%.M  alpha   hy   hz   Ty+  Ty-  Tz+  Tz-
1   1e-4    1.5  1.5  80   20   30  -10
0                               % number of joints with support settlements
				% End   Static Load Case 2 of 3

				% Begin Static Load Case 3 of 3
0				% number of loaded joints
0                               % number of distributed loads
2                               % number of internal concentrated loads
%.M    Px   Py   Pz   x    
  1    0    10   -90  30
  2    0   -20    20  30
0                               % number of members with temperature loads
0                               % number of joints with support settlements
				% End   Static Load Case 3 of 3


6				% number of desired dynamic modes of vibration
1                               % 1: subspace Jacobi     2: Stodola
0				% 0: consistent mass ... 1: lumped mass matrix
1e-4				% mode shape tolerance
0.0				% shift value ... for unrestrained structures

1       7.e-5   0.0		% beam numbers, density, and extra beam mass 
2       7.e-5   0.0
3       7.e-5   0.0
4       7.e-5   0.0

0                               % number of joints with extra inertia
%.j      M      Ixx     Iyy     Izz -  joints and concentrated mass and inertia

6				% number of modes to animate
 1  2  3  4 5 6 		% modes to animate
1                               % pan during animation


