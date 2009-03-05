Example C: a tetrahedral frame - static and dynamic analysis 

18                      % number of joints 

%.joint  x       y       z       r
%       in      in      in      in

1	0.0	0.0	0.0	0.0
2	0.0	100.0	0.0	0.0
3	0.0	50.0	100.0	0.0
4	100.0	0.0	0.0	0.0
5	100.0	100.0	0.0	0.0
6	100.0	50.0	100.0	0.0
7	200.0	0.0	0.0	0.0
8	200.0	100.0	0.0	0.0
9	200.0	50.0	100.0	0.0
10	300.0	0.0	0.0	0.0
11	300.0	100.0	0.0	0.0
12	300.0	50.0	100.0	0.0
13	400.0	0.0	0.0	0.0
14	400.0	100.0	0.0	0.0
15	400.0	50.0	100.0	0.0
16	500.0	0.0	0.0	0.0
17	500.0	100.0	0.0	0.0
18	500.0	50.0	100.0	0.0

4                               % number of joints with reactions
%.J     x y z xx yy zz          1=fixed, 0=free
 1	1 1 1  0  0  0
 3	1 1 1  0  0  0
 16	1 1 1  0  0  0
 18	1 1 1  0  0  0

48 			% number of members
%.m j1 j2        Ax     Asy     Asz     Jxx     Iyy     Izz     E       G    p
%               in^2    in^2    in^2    in^4    in^4    in^4    ksi     ksi  deg

1   1  2	10.0	8.0	6.0	50.0	100.0	200.0	29000  	11500 0
2   2  3	10.0	7.0	5.0	50.0	100.0	200.0	29000  	11500 0
3   1  3	20.0	8.0	5.0	50.0	100.0	200.0	29000  	11500 0
4   1  4	10.0	8.0	6.0	50.0	100.0	200.0	29000  	11500 0
5   2  5	10.0	7.0	5.0	50.0	100.0	200.0	29000   11500 0
6   3  6	20.0	8.0	5.0	50.0	100.0	200.0	29000  	11500 0
7   3  4	10.0	8.0	6.0	50.0	100.0	200.0	29000  	11500 0
8   2  4	10.0	7.0	5.0	50.0	100.0	200.0	29000  	11500 0
9   2  6	20.0	8.0	5.0	50.0	100.0	200.0	29000  	11500 0
10  4  5	10.0	8.0	6.0	50.0	100.0	200.0	29000  	11500 0
11  5  6	10.0	7.0	5.0	50.0	100.0	200.0	29000  	11500 0
12  4  6	20.0	8.0	5.0	50.0	100.0	200.0	29000  	11500 0
13  4  7	10.0	8.0	6.0	50.0	100.0	200.0	29000  	11500 0
14  5  8	10.0	7.0	5.0	50.0	100.0	200.0	29000  	11500 0
15  6  9	20.0	8.0	5.0	50.0	100.0	200.0	29000  	11500 0
16  4  9	10.0	8.0	6.0	50.0	100.0	200.0	29000	11500 0
17  5  9	10.0	7.0	5.0	50.0	100.0	200.0	29000  	11500 0
18  4  8	20.0	8.0	5.0	50.0	100.0	200.0	29000  	11500 0
19  7  8	10.0	8.0	6.0	50.0	100.0	200.0	29000  	11500 0
20  8  9	10.0	7.0	5.0	50.0	100.0	200.0	29000  	11500 0
21  7  9	20.0	8.0	5.0	50.0	100.0	200.0	29000  	11500 0
22  7  10	10.0	8.0	6.0	50.0	100.0	200.0	29000  	11500 0
23  8  11	10.0	7.0	5.0	50.0	100.0	200.0	29000  	11500 0
24  9  12	20.0	8.0	5.0	50.0	100.0	200.0	29000  	11500 0
25  9  10	10.0	8.0	6.0	50.0	100.0	200.0	29000  	11500 0
26  8  10	10.0	7.0	5.0	50.0	100.0	200.0	29000  	11500 0
27  8  12	20.0	8.0	5.0	50.0	100.0	200.0	29000  	11500 0
28  10 11	10.0	8.0	6.0	50.0	100.0	200.0	29000  	11500 0
29  11 12	10.0	7.0	5.0	50.0	100.0	200.0	29000  	11500 0
30  10 12	20.0	8.0	5.0	50.0	100.0	200.0	29000  	11500 0
31  10 13	10.0	8.0	6.0	50.0	100.0	200.0	29000  	11500 0
32  11 14	10.0	7.0	5.0	50.0	100.0	200.0	29000  	11500 0
33  12 15	20.0	8.0	5.0	50.0	100.0	200.0	29000  	11500 0
34  10 15	10.0	8.0	6.0	50.0	100.0	200.0	29000  	11500 0
35  11 15	10.0	7.0	5.0	50.0	100.0	200.0	29000  	11500 0
36  10 14	20.0	8.0	5.0	50.0	100.0	200.0	29000  	11500 0
37  13 14	10.0	8.0	6.0	50.0	100.0	200.0	29000  	11500 0
38  14 15	10.0	7.0	5.0	50.0	100.0	200.0	29000  	11500 0
39  13 15	20.0	8.0	5.0	50.0	100.0	200.0	29000  	11500 0
40  13 16	10.0	8.0	6.0	50.0	100.0	200.0	29000  	11500 0
41  14 17	10.0	7.0	5.0	50.0	100.0	200.0	29000  	11500 0
42  15 18	20.0	8.0	5.0	50.0	100.0	200.0	29000  	11500 0
43  15 16	10.0	8.0	6.0	50.0	100.0	200.0	29000  	11500 0
44  14 16	10.0	7.0	5.0	50.0	100.0	200.0	29000  	11500 0
45  14 18	20.0	8.0	5.0	50.0	100.0	200.0	29000  	11500 0
46  16 17	10.0	8.0	6.0	50.0	100.0	200.0	29000	11500 0
47  17 18	10.0	7.0	5.0	50.0	100.0	200.0	29000	11500 0
48  16 18	20.0	8.0	5.0	50.0	100.0	200.0	29000	11500 0

1                               % 1: include shear deformation
1                               % 1: include geometric stiffness
2.0                             % exaggerate mesh deformations


 1                      % number of static load cases
				% Begin Static Load Case 1 of 1
0                               % number of loaded joints

5				% number of uniform loads
%.j      Ux       Uy     Uz
%       k/in     k/in   k/in
  5	0.0	-10.0	0.0
 14	0.0	-10.0	0.0
 23	0.0	-10.0	0.0
 32	0.0	-10.0	0.0
 41	0.0	-10.0	0.0

0                               % number of trapezoidal loads
0                               % number of internal concentrated loads
0                               % number of temperature loads
0                               % number of support settlements
				% End  Static Load Case 1 of 1

10                              % number of desired dynamic modes of vibration
1                               % 1: subspace Jacobi     2: Stodola
0                               % 0: consistent mass ... 1: lumped mass matrix
1e-9                            % mode shape tolerance
0.0                             % shift value ... for unrestrained structures


%.M     density mass
%       k/in^3  kip
1	7.32e-7	0.0	        % bar numbers, density, and extra mass
2	7.32e-7	0.0
3	7.32e-7	0.0
4	7.32e-7	0.0
5	7.32e-7	0.0
6	7.32e-7	0.0
7	7.32e-7	0.0
8	7.32e-7	0.0
9	7.32e-7	0.0
10	7.32e-7	0.0
11	7.32e-7	0.0
12	7.32e-7	0.0
13	7.32e-7	0.0
14	7.32e-7	0.0
15	7.32e-7	0.0
16	7.32e-7	0.0
17	7.32e-7	0.0
18	7.32e-7	0.0
19	7.32e-7	0.0
20	7.32e-7	0.0
21	7.32e-7	0.0
22	7.32e-7	0.0
23	7.32e-7	0.0
24	7.32e-7	0.0
25	7.32e-7	0.0
26	7.32e-7	0.0
27	7.32e-7	0.0
28	7.32e-7	0.0
29	7.32e-7	0.0
30	7.32e-7	0.0
31	7.32e-7	0.0
32	7.32e-7	0.0
33	7.32e-7	0.0
34	7.32e-7	0.0
35	7.32e-7	0.0
36	7.32e-7	0.0
37	7.32e-7	0.0
38	7.32e-7	0.0
39	7.32e-7	0.0
40	7.32e-7	0.0
41	7.32e-7	0.0
42	7.32e-7	0.0
43	7.32e-7	0.0
44	7.32e-7	0.0
45	7.32e-7	0.0
46	7.32e-7	0.0
47	7.32e-7	0.0
48	7.32e-7	0.0

0                         % number of joints with extra mass or inertia
%.j      M      Ixx     Iyy     Izz -  joints and concentrated mass and inertia

5                               % number of modes to animate
 1  2  3  4   5                 % modes to animate
2                               % pan rate during animation

