Example J: tesseract 

16			% number of joints
% joint    x         y         z         r
 
 1       100       100      -100         0
 2      -100       100      -100         0
 3      -100      -100      -100         0         
 4       100      -100      -100         0

 5       100       100       100         0
 6      -100       100       100         0
 7      -100      -100       100         0         
 8       100      -100       100         0

 9        70        70       -70         0
10       -70        70       -70         0
11       -70       -70       -70         0         
12        70       -70       -70         0

13        70        70        70         0
14       -70        70        70         0
15       -70       -70        70         0         
16        70       -70        70         0

0			% number of joints with reactions

32			% number of members
% m    j1 j2   Ax  Asy   Asz   Jxx    Iyy   Izz       E      G p

 1      1  2  100  60    60    500   1000  1000   10000   7000 0
 2      2  3  100  60    60    500   1000  1000   10000   7000 0
 3      3  4  100  60    60    500   1000  1000   10000   7000 0
 4      4  1  100  60    60    500   1000  1000   10000   7000 0

 5      5  6  100  60    60    500   1000  1000   10000   7000 0
 6      6  7  100  60    60    500   1000  1000   10000   7000 0
 7      7  8  100  60    60    500   1000  1000   10000   7000 0
 8      8  5  100  60    60    500   1000  1000   10000   7000 0

 9      9 10  100  60    60    500   1000  1000   10000   7000 0
10     10 11  100  60    60    500   1000  1000   10000   7000 0
11     11 12  100  60    60    500   1000  1000   10000   7000 0
12     12  9  100  60    60    500   1000  1000   10000   7000 0

13     13 14  100  60    60    500   1000  1000   10000   7000 0
14     14 15  100  60    60    500   1000  1000   10000   7000 0
15     15 16  100  60    60    500   1000  1000   10000   7000 0
16     16 13  100  60    60    500   1000  1000   10000   7000 0

17      1  5  100  60    60    500   1000  1000   10000   7000 0
18      2  6  100  60    60    500   1000  1000   10000   7000 0
19      3  7  100  60    60    500   1000  1000   10000   7000 0
20      4  8  100  60    60    500   1000  1000   10000   7000 0

21      9 13  100  60    60    500   1000  1000   10000   7000 0
22     10 14  100  60    60    500   1000  1000   10000   7000 0
23     11 15  100  60    60    500   1000  1000   10000   7000 0
24     12 16  100  60    60    500   1000  1000   10000   7000 0

25      1  9  100  60    60    500   1000  1000   10000   7000 0
26      2 10  100  60    60    500   1000  1000   10000   7000 0
27      3 11  100  60    60    500   1000  1000   10000   7000 0
28      4 12  100  60    60    500   1000  1000   10000   7000 0

29      5 13  100  60    60    500   1000  1000   10000   7000 0
30      6 14  100  60    60    500   1000  1000   10000   7000 0
31      7 15  100  60    60    500   1000  1000   10000   7000 0
32      8 16  100  60    60    500   1000  1000   10000   7000 0

1                               % 1: include shear deformation
1                               % 1: include geometric stiffness
5.0                             % exaggerate mesh deformations
1                               % 1: stiffness analysis, 0: data check only

1			% number of static load cases
				% Begin Static Load Case 1 of 1
0                               % number of loaded joints
0                               % number of distributed loads
0                               % number of internal concentrated loads
0                               % number of members with temperature loads
0                               % number of joints with support settlements
				% End   Static Load Case 1 of 1

22                              % number of desired dynamic modes of vibration
1                               % 1: subspace Jacobi     2: Stodola
0                               % 0: consistent mass ... 1: lumped mass matrix
1e-6                            % mode shape tolerance
1.0                             % shift value ... for unrestrained structures

 1     3e-7  0			% bar numbers, density, and extra mass
 2     3e-7  0
 3     3e-7  0
 4     3e-7  0
 5     3e-7  0
 6     3e-7  0
 7     3e-7  0
 8     3e-7  0
 9     3e-7  0
10     3e-7  0
11     3e-7  0
12     3e-7  0
13     3e-7  0
14     3e-7  0
15     3e-7  0
16     3e-7  0
17     3e-7  0
18     3e-7  0
19     3e-7  0
20     3e-7  0
21     3e-7  0
22     3e-7  0
23     3e-7  0
24     3e-7  0
25     3e-7  0
26     3e-7  0
27     3e-7  0
28     3e-7  0
29     3e-7  0
30     3e-7  0
31     3e-7  0
32     3e-7  0

0                % number of joints with extra mass or inertia

8                               % number of modes to animate
 7  10  12  13 15  18  20  22   % modes to animate
1                               % pan during animation


