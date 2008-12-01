pyramid-shaped frame --- static and dynamic analysis

5			% number of joints 
4			% number of members			
3			% number of static load cases

% joint  x       y       z       r

1	0.0	0.0	100.0	0.0
2	-120.0	-90.0	0.0	0.0	
3	 120.0	-90.0	0.0	0.0	
4	 120.0	 90.0	0.0	0.0	
5	-120.0	 90.0	0.0	0.0	

% m j1 j2 Ax    Asy     Asz     Jxx     Iyy     Izz         E      G   p

1 1 2	10.	8.0	8.0	500.0	300.0	200.0	1000.0	700.0  0
2 1 3	10.	8.0	8.0	500.0	300.0	200.0	1000.0	700.0  0
3 1 4	10.	8.0	8.0	500.0	300.0	200.0	1000.0	700.0  0
4 1 5	10.	8.0	8.0	500.0	300.0	200.0	1000.0	700.0  0

1                               % 1: include shear deformation
0                               % 1: include geometric stiffness
/tmp/exB-msh                    % mesh data file name
exB.plt                         % plot file name
10.0                            % exaggerate mesh deformations
1                               % 1: stiffness analysis, 0: data check only

				% Begin Static Load Case 1 of 3
1				% number of loaded joints
%J      Fx       Fy     Fz      Mxx     Myy     Mzz
 1	10.00	-20.00	-100	0.0	0.0	0.0
0                               % number of distributed loads
0                               % number of internal concentrated loads
0                               % number of members with temperature loads
				% End   Static Load Case 1 of 3

				% Begin Static Load Case 2 of 3
0				% number of loaded joints
2                               % number of distributed loads
% M    Wx   Wy   Wz
  2    0    0.1  0
  1    0    0    0.1
0                               % number of internal concentrated loads
1                               % number of members with temperature loads
%M  alpha   hy   hz   Ty+  Ty-  Tz+  Tz-
1   1e-4    1.5  1.5  80   20   30  -10
				% End   Static Load Case 2 of 3

				% Begin Static Load Case 3 of 3
0				% number of loaded joints
0                               % number of distributed loads
2                               % number of internal concentrated loads
% M    Px   Py   Pz   x    
  1    0    10   -90  30
  2    0   -20    20  30
0                               % number of members with temperature loads
				% End   Static Load Case 3 of 3

4                               % number of joints with reactions
% J     x y z xx yy zz          1= fixed, 0=free

  2	1 1 1 1 1 1
  3	1 1 1 1 1 1
  4	1 1 1 1 1 1
  5	1 1 1 1 1 1
 
0                               % number of joints with support settlements

6				% number of desired dynamic modes of vibration
1                               % 1: subspace Jacobi     2: Stodola
0				% 0: consistent mass ... 1: lumped mass matrix
/tmp/exB-m			% mode shape data file
1e-4				% mode shape tolerance
0.0				% shift value ... for unrestrained structures

1       7.e-5   0.0		% beam numbers, density, and extra beam mass 
2       7.e-5   0.0
3       7.e-5   0.0
4       7.e-5   0.0

0                               % number of joints with extra inertia
% j      M      Ixx     Iyy     Izz -  joints and concentrated mass and inertia

6				% number of modes to animate
 1  2  3  4 5 6 		% modes to animate
1                               % pan during animation



________________________________________________________________________________
FRAME3DD version: 20080909               http://frame3dd.sf.net/
GPL Copyright (C) 1992-2008, Henri P. Gavin 
FRAME3DD is distributed in the hope that it will be useful but with no warranty.
For details see the GNU Public Licence: http://www.fsf.org/copyleft/gpl.html
________________________________________________________________________________

pyramid-shaped frame --- static and dynamic analysis 
Mon Dec  1 13:06:28 2008
________________________________________________________________________________
5 JOINTS;    4 MEMBERS;    3 LOAD CASES;
4 FIXED JOINTS;   0 PRESCRIBED DISPLACEMENTS;
For 2D problems, the Y-axis is vertical. 
For 3D problems, the Z-axis is vertical. 
________________________________________________________________________________
J O I N T   D A T A                                         R E S T R A I N T S
  Joint      X              Y              Z         radius  Fx Fy Fz Mx My Mz
    1       0.000000       0.000000     100.000000    0.000   0  0  0  0  0  0
    2    -120.000000     -90.000000       0.000000    0.000   1  1  1  1  1  1
    3     120.000000     -90.000000       0.000000    0.000   1  1  1  1  1  1
    4     120.000000      90.000000       0.000000    0.000   1  1  1  1  1  1
    5    -120.000000      90.000000       0.000000    0.000   1  1  1  1  1  1
M E M B E R   D A T A							(local)
  Member J1    J2     Ax   Asy   Asz    Jxx     Iyy     Izz       E       G roll
    1     1     2   10.0   8.0   8.0  500.0   300.0   200.0   1000.0   700.0   0
    2     1     3   10.0   8.0   8.0  500.0   300.0   200.0   1000.0   700.0   0
    3     1     4   10.0   8.0   8.0  500.0   300.0   200.0   1000.0   700.0   0
    4     1     5   10.0   8.0   8.0  500.0   300.0   200.0   1000.0   700.0   0
  Include shear deformations.
  Neglect geometric stiffness.

L O A D   C A S E   1   O F   3  ... 

 1 joints  with concentrated loads
 0 members with uniformly distributed loads
 0 members with concentrated point loads
 0 members with temperature loads
 J O I N T   L O A D S  +  E Q U I V A L E N T   J O I N T   L O A D S  (global)
  Joint       Fx          Fy          Fz          Mxx         Myy         Mzz
     1      10.000     -20.000    -100.000       0.000       0.000       0.000

L O A D   C A S E   2   O F   3  ... 

 0 joints  with concentrated loads
 2 members with uniformly distributed loads
 0 members with concentrated point loads
 1 members with temperature loads
 J O I N T   L O A D S  +  E Q U I V A L E N T   J O I N T   L O A D S  (global)
  Joint       Fx          Fy          Fz          Mxx         Myy         Mzz
     1      21.378      19.188      24.141    -877.323     500.272     890.987
     2     -23.969     -17.977      -9.141     997.508    -590.411    -665.640
     3       5.408       7.211       0.000    -120.185      90.139    -225.347
 U N I F O R M   M E M B E R   L O A D S					(local)
  Member      Wx               Wy               Wz
     2       0.00000000       0.10000000       0.00000000
     1       0.00000000       0.00000000       0.10000000
 M E M B E R   T E M P E R A T U R E   C H A N G E S			(local)
  Member    coef      hy        hz        Ty+       Ty-       Tz+       Tz-
     1  1.00e-04     1.500     1.500    80.000    20.000    30.000   -10.000

L O A D   C A S E   3   O F   3  ... 

 0 joints  with concentrated loads
 0 members with uniformly distributed loads
 2 members with concentrated point loads
 0 members with temperature loads
 J O I N T   L O A D S  +  E Q U I V A L E N T   J O I N T   L O A D S  (global)
  Joint       Fx          Fy          Fz          Mxx         Myy         Mzz
     1      35.737      -0.601     -48.621     361.106   -1067.740    -105.385
     2       7.586       3.632     -12.371    -614.576     933.536    -102.692
     3      -0.509      -3.734       2.749     258.281     116.171     205.384
 C O N C E T R A T E D   P O I N T   L O A D S				(local)
  Member      Px          Py          Pz          x
     1       0.000      10.000     -90.000      30.000
     2       0.000     -20.000      20.000      30.000

E L A S T I C   S T I F F N E S S   A N A L Y S I S   via  L D L'  decomposition


L O A D   C A S E   1   O F   3  ... 

J O I N T   D I S P L A C E M E N T S					(global)
  Joint    X-dsp       Y-dsp       Z-dsp       X-rot       Y-rot       Z-rot
     1    0.101193   -0.356839   -1.429811    0.002076    0.000520    0.0     
M E M B E R   E N D   F O R C E S					(local)
  Member Joint      Nx          Vy         Vz         Txx        Myy        Mzz
     1      1     50.140c      0.100     -0.717     -3.187     66.044      7.806
     1      2    -50.140c     -0.100      0.717      3.187     63.282     10.235
     2      1     57.612c     -0.064     -0.708      2.180     66.624     -4.902
     2      3    -57.612c      0.064      0.708     -2.180     61.093     -6.562
     3      1     37.849c     -0.100     -0.716      3.187     63.179     -7.806
     3      4    -37.849c      0.100      0.716     -3.187     65.941    -10.235
     4      1     30.376c      0.064     -0.725     -2.180     62.600      4.902
     4      5    -30.376c     -0.064      0.725      2.180     68.130      6.562
R E A C T I O N S							(global)
  Joint       Fx          Fy          Fz         Mxx         Myy         Mzz
     2      32.996      24.872      28.409      31.306     -55.623       6.748
     3     -37.996      28.577      32.547      32.293      52.146      -4.251
     4     -24.936     -18.577      21.591     -46.228      47.755      -6.748
     5      19.936     -14.872      17.453     -45.241     -51.232       4.251
R M S   E Q U I L I B R I U M    E R R O R: 8.011e-17

L O A D   C A S E   2   O F   3  ... 

J O I N T   D I S P L A C E M E N T S					(global)
  Joint    X-dsp       Y-dsp       Z-dsp       X-rot       Y-rot       Z-rot
     1    0.242338    0.433439    0.345171   -0.056188    0.025442    0.061201
M E M B E R   E N D   F O R C E S					(local)
  Member Joint      Nx          Vy         Vz         Txx        Myy        Mzz
     1      1     -1.571t      2.374     -6.056    -17.957    714.258   -511.211
     1      2      1.571t     -2.374    -11.972     17.957  -1247.464    939.266
     2      1    -13.676t     -8.176      0.876   -163.181   -101.233   -175.905
     2      3     13.676t     -9.851     -0.876    163.181    -56.772    326.886
     3      1     10.330c      1.337     -2.612   -113.862    325.379    158.749
     3      4    -10.330c     -1.337      2.612    113.862    145.435     82.297
     4      1     -7.566t      2.874     -0.530     31.362     70.037    352.609
     4      5      7.566t     -2.874      0.530    -31.362     25.576    165.511
R E A C T I O N S							(global)
  Joint       Fx          Fy          Fz         Mxx         Myy         Mzz
     2       2.842       5.099     -10.833   -1177.240     676.400     771.556
     3       2.803     -14.417      -8.315     219.615    -235.676     181.469
     4      -4.915      -5.358       7.903      25.050     200.582       5.316
     5      -3.547       6.253      -3.755     -67.917      18.968     155.110
R M S   E Q U I L I B R I U M    E R R O R: 1.201e-16

L O A D   C A S E   3   O F   3  ... 

J O I N T   D I S P L A C E M E N T S					(global)
  Joint    X-dsp       Y-dsp       Z-dsp       X-rot       Y-rot       Z-rot
     1    0.306117   -0.048209   -0.695194    0.022369   -0.050071   -0.007239
M E M B E R   E N D   F O R C E S					(local)
  Member Joint      Nx          Vy         Vz         Txx        Myy        Mzz
     1      1     11.423c     -8.238     71.806     27.418   -746.401   -115.416
     1      2    -11.423c     -1.762     18.194    -27.418   1326.305    133.075
     2      1     34.028c     17.516    -15.506     85.234    100.698    348.984
     2      3    -34.028c      2.484     -4.494    -85.234   -310.808   -196.807
     3      1     31.358c     -0.555      2.628    -11.827   -325.887    -64.175
     3      4    -31.358c      0.555     -2.628     11.827   -147.900    -35.861
     4      1      8.753c     -1.247     -1.887    -69.643    214.381   -148.605
     4      5     -8.753c      1.247      1.887     69.643    125.733    -76.240
R E A C T I O N S							(global)
  Joint       Fx          Fy          Fz         Mxx         Myy         Mzz
     2      -1.527       1.057      21.474     754.980   -1091.646     125.934
     3     -23.154      20.471      15.137    -330.555    -140.594    -116.474
     4     -22.373     -16.086      15.208      80.699    -124.351     -36.398
     5       4.241      -4.740       6.425     -87.965     -91.193    -102.066
R M S   E Q U I L I B R I U M    E R R O R: 1.948e-16

M O D A L   A N A L Y S I S   R E S U L T S
  Total Mass:  5.047772e-01     Structural Mass:  5.047772e-01 
J O I N T   M A S S E S	(diagonal of the mass matrix)			(global)
  Joint X-mass      Y-mass      Z-mass      X-inrta     Y-inrta     Z-inrta
     1 1.79213e-01 1.82997e-01 1.81959e-01 9.17335e+01 1.20839e+02 1.11687e+02
     2 3.24803e+02 3.24803e+02 3.24803e+02 3.24803e+02 3.24803e+02 3.24803e+02
     3 3.24803e+02 3.24803e+02 3.24803e+02 3.24803e+02 3.24803e+02 3.24803e+02
     4 3.24803e+02 3.24803e+02 3.24803e+02 3.24803e+02 3.24803e+02 3.24803e+02
     5 3.24803e+02 3.24803e+02 3.24803e+02 3.24803e+02 3.24803e+02 3.24803e+02
  Use consistent mass matrix.
N A T U R A L   F R E Q U E N C I E S   & 
M A S S   N O R M A L I Z E D   M O D E   S H A P E S 
 convergence tolerance: 1.000e-04 
  MODE     1:   f= 1.816585 Hz,  T= 0.550484 sec
		X- modal participation factor =   1.3449e-08 
		Y- modal participation factor =  -1.7942e-08 
		Z- modal participation factor =   1.6169e-10 
  Joint   X-dsp       Y-dsp       Z-dsp       X-rot       Y-rot       Z-rot
     1  -5.996e-09   5.963e-09  -6.109e-11   2.468e-10   1.021e-10   9.457e-02
     2   8.412e-07  -1.122e-06   2.470e-13   2.694e-05   2.021e-05  -4.799e-05
     3   8.412e-07   1.122e-06  -8.797e-13  -2.694e-05   2.020e-05  -4.799e-05
     4  -8.412e-07   1.122e-06   2.785e-13  -2.694e-05  -2.020e-05  -4.799e-05
     5  -8.412e-07  -1.122e-06   8.659e-13   2.694e-05  -2.020e-05  -4.799e-05
  MODE     2:   f= 1.927166 Hz,  T= 0.518897 sec
		X- modal participation factor =  -4.8240e-11 
		Y- modal participation factor =   5.2790e-01 
		Z- modal participation factor =  -1.4977e-08 
  Joint   X-dsp       Y-dsp       Z-dsp       X-rot       Y-rot       Z-rot
     1   1.097e-08   7.745e-01   9.956e-09   7.860e-02   5.544e-10   6.081e-10
     2   3.546e-08   1.259e-06  -7.566e-07  -4.435e-05   2.275e-05   3.558e-05
     3  -3.545e-08   1.259e-06  -7.566e-07  -4.434e-05  -2.274e-05  -3.557e-05
     4   3.545e-08   1.259e-06   7.566e-07  -4.434e-05   2.274e-05  -3.557e-05
     5  -3.546e-08   1.259e-06   7.566e-07  -4.434e-05  -2.274e-05   3.557e-05
  MODE     3:   f= 2.008918 Hz,  T= 0.497780 sec
		X- modal participation factor =  -4.3969e-01 
		Y- modal participation factor =   1.3859e-09 
		Z- modal participation factor =   2.8423e-08 
  Joint   X-dsp       Y-dsp       Z-dsp       X-rot       Y-rot       Z-rot
     1  -3.665e-01  -1.840e-08  -1.690e-08   1.112e-09   8.214e-02   5.006e-10
     2  -1.202e-06  -1.798e-08   1.172e-06   2.585e-05  -6.011e-05   2.550e-05
     3  -1.202e-06   1.798e-08  -1.171e-06  -2.585e-05  -6.011e-05   2.549e-05
     4  -1.202e-06  -1.798e-08  -1.171e-06   2.585e-05  -6.011e-05  -2.549e-05
     5  -1.202e-06   1.798e-08   1.171e-06  -2.585e-05  -6.011e-05  -2.549e-05
  MODE     4:   f= 3.119394 Hz,  T= 0.320575 sec
		X- modal participation factor =   8.5668e-08 
		Y- modal participation factor =   7.1311e-08 
		Z- modal participation factor =   5.9565e-01 
  Joint   X-dsp       Y-dsp       Z-dsp       X-rot       Y-rot       Z-rot
     1  -1.043e-07  -1.580e-07   2.343e+00   1.517e-08  -1.295e-08  -6.892e-11
     2   3.252e-07   2.439e-07   3.134e-06   6.239e-05  -8.319e-05  -1.465e-10
     3  -3.250e-07   2.438e-07   3.133e-06   6.236e-05   8.315e-05  -1.166e-09
     4  -3.250e-07  -2.437e-07   3.133e-06  -6.236e-05   8.315e-05  -1.601e-10
     5   3.250e-07  -2.437e-07   3.133e-06  -6.236e-05  -8.314e-05   1.184e-09
  MODE     5:   f= 4.021047 Hz,  T= 0.248691 sec
		X- modal participation factor =   4.3050e-07 
		Y- modal participation factor =   2.7858e-01 
		Z- modal participation factor =   2.5935e-06 
  Joint   X-dsp       Y-dsp       Z-dsp       X-rot       Y-rot       Z-rot
     1  -2.868e-07   2.962e+00  -7.805e-07  -1.118e-01  -5.346e-08   1.296e-08
     2   6.190e-07   1.005e-06   5.548e-06   8.010e-05  -1.454e-04   1.665e-05
     3  -6.180e-07   1.003e-06   5.540e-06   8.000e-05   1.452e-04  -1.661e-05
     4   6.182e-07   1.003e-06  -5.540e-06   8.000e-05  -1.452e-04  -1.661e-05
     5  -6.181e-07   1.003e-06  -5.540e-06   8.001e-05   1.452e-04   1.661e-05
  MODE     6:   f= 4.706763 Hz,  T= 0.212460 sec
		X- modal participation factor =   4.1314e-01 
		Y- modal participation factor =   1.0282e-06 
		Z- modal participation factor =   6.3788e-06 
  Joint   X-dsp       Y-dsp       Z-dsp       X-rot       Y-rot       Z-rot
     1   2.848e+00  -5.367e-07  -1.347e-06   1.264e-07   7.391e-02  -5.890e-08
     2   4.162e-06   8.268e-07   7.098e-06   1.339e-04  -8.991e-05  -6.721e-05
     3   4.149e-06  -8.247e-07  -7.078e-06  -1.335e-04  -8.971e-05  -6.699e-05
     4   4.149e-06   8.249e-07  -7.078e-06   1.335e-04  -8.968e-05   6.699e-05
     5   4.149e-06  -8.243e-07   7.078e-06  -1.335e-04  -8.969e-05   6.700e-05
M A T R I X    I T E R A T I O N S: 2
There are 6 modes below 4.706764 Hz. ... All 6 modes were found.
