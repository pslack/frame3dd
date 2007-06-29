pyramid-shaped frame --- static and dynamic analysis

5	4			% number of joints and number of members

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

1				% number of loaded joints
%J      Fx       Fy     Fz      Mxx     Myy     Mzz
 1	10.00	-20.00	-280	0.0	0.0	0.0

0                               % number of distributed loads
0                               % number of internal concentrated loads
0                               % number of members with temperature loads


4                               % number of joints with reactions
% J     x y z xx yy zz          1= fixed, 0=free

  2	1 1 1 1 1 1
  3	1 1 1 1 1 1
  4	1 1 1 1 1 1
  5	1 1 1 1 1 1
 
0                               % number of joints with support settlements

6				% number of desired dynamic modes of vibration
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
 1  2  3  4 5 6 			% modes to animate
1                               % pan during animation



________________________________________________________________________________
-- FRAME version:   1 Mar 2007, GPL Copyright (C) 1992-2007, Henri P. Gavin --
                     http://www.duke.edu/~hpgavin/frame/ 
 FRAME is distributed in the hope that it will be useful but with no warranty;
 for details see the GNU Public Licence: http://www.fsf.org/copyleft/gpl.html
________________________________________________________________________________

pyramid-shaped frame --- static and dynamic analysis 
Thu Mar  1 17:31:04 2007
________________________________________________________________________________
JOINTS: 5    MEMBERS: 4   FIXED JOINTS: 4   PRESCRIBED DISPLACEMENTS: 0
JOINT LOADS: 1   UNIFORM MEMBER LOADS: 0   CONCENTRATED MEMBER LOADS: 0   

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
J O I N T   L O A D S  +  E Q U I V A L E N T   J O I N T   L O A D S	(global)
  Joint       Fx          Fy          Fz          Mxx         Myy         Mzz
     1      10.000     -20.000    -280.000       0.000       0.000       0.000

E L A S T I C   S T I F F N E S S   A N A L Y S I S   via  L D L'  decomposition

J O I N T   D I S P L A C E M E N T S					(global)
  Joint    X-dsp       Y-dsp       Z-dsp       X-rot       Y-rot       Z-rot
     1    0.101193   -0.356839   -4.003471    0.002076    0.000520    0.0     
M E M B E R   E N D   F O R C E S					(local)
  Member Joint      Nx          Vy         Vz         Txx        Myy        Mzz
     1      1    129.329c      0.100     -2.008     -3.187    182.345      7.806
     1      2   -129.329c     -0.100      2.008      3.187    179.583     10.235
     2      1    136.802c     -0.064     -1.999      2.180    182.925     -4.902
     2      3   -136.802c      0.064      1.999     -2.180    177.394     -6.562
     3      1    117.038c     -0.100     -2.006      3.187    179.480     -7.806
     3      4   -117.038c      0.100      2.006     -3.187    182.242    -10.235
     4      1    109.566c      0.064     -2.015     -2.180    178.900      4.902
     4      5   -109.566c     -0.064      2.015      2.180    184.431      6.562
R E A C T I O N S							(global)
  Joint       Fx          Fy          Fz         Mxx         Myy         Mzz
     2      85.136      63.977      73.409     101.086    -148.664       6.748
     3     -90.136      67.681      77.547     102.073     145.187      -4.251
     4     -77.075     -57.681      66.591    -116.009     140.796      -6.748
     5      72.075     -53.977      62.453    -115.022    -144.273       4.251
R M S   E Q U I L I B R I U M    E R R O R: 1.264e-16

M O D A L   A N A L Y S I S   R E S U L T S
  Total Mass:  5.047772e-01     Structural Mass:  5.047772e-01 
J O I N T   M A S S E S	(diagonal of the mass matrix)			(global)
  Joint X-mass      Y-mass      Z-mass      X-inrta     Y-inrta     Z-inrta
     1 1.79213e-01 1.82997e-01 1.81959e-01 9.17335e+01 1.20839e+02 1.11687e+02
     2 6.49607e+02 6.49607e+02 6.49607e+02 6.49607e+02 6.49607e+02 6.49607e+02
     3 6.49607e+02 6.49607e+02 6.49607e+02 6.49607e+02 6.49607e+02 6.49607e+02
     4 6.49607e+02 6.49607e+02 6.49607e+02 6.49607e+02 6.49607e+02 6.49607e+02
     5 6.49607e+02 6.49607e+02 6.49607e+02 6.49607e+02 6.49607e+02 6.49607e+02
  Use consistent mass matrix.
N A T U R A L   F R E Q U E N C I E S   & 
M A S S   N O R M A L I Z E D   M O D E   S H A P E S 
 convergence tolerance: 1.000e-04 
  MODE     1:   f= 1.816581 Hz,  T= 0.550485 sec
		X- modal participation factor =   1.2629e-07 
		Y- modal participation factor =  -1.6865e-07 
		Z- modal participation factor =   7.6530e-10 
  Joint   X-dsp       Y-dsp       Z-dsp       X-rot       Y-rot       Z-rot
     1  -2.732e-08   2.657e-08  -1.439e-10   1.200e-09   5.087e-10   9.457e-02
     2   8.478e-07  -1.130e-06   5.716e-13   2.715e-05   2.036e-05  -4.837e-05
     3   8.476e-07   1.130e-06  -1.951e-12  -2.715e-05   2.036e-05  -4.836e-05
     4  -8.476e-07   1.130e-06   6.399e-13  -2.715e-05  -2.036e-05  -4.836e-05
     5  -8.476e-07  -1.130e-06   1.935e-12   2.714e-05  -2.036e-05  -4.836e-05
  MODE     2:   f= 1.927161 Hz,  T= 0.518898 sec
		X- modal participation factor =   3.4103e-09 
		Y- modal participation factor =   5.2956e-01 
		Z- modal participation factor =  -1.4294e-07 
  Joint   X-dsp       Y-dsp       Z-dsp       X-rot       Y-rot       Z-rot
     1   5.112e-08   7.745e-01   4.604e-08   7.860e-02   2.680e-09   2.806e-09
     2   3.577e-08   1.270e-06  -7.633e-07  -4.474e-05   2.295e-05   3.589e-05
     3  -3.576e-08   1.269e-06  -7.631e-07  -4.472e-05  -2.294e-05  -3.588e-05
     4   3.576e-08   1.269e-06   7.631e-07  -4.472e-05   2.294e-05  -3.588e-05
     5  -3.576e-08   1.269e-06   7.631e-07  -4.472e-05  -2.294e-05   3.588e-05
  MODE     3:   f= 2.008911 Hz,  T= 0.497782 sec
		X- modal participation factor =  -4.4128e-01 
		Y- modal participation factor =   4.1944e-09 
		Z- modal participation factor =   2.6587e-07 
  Joint   X-dsp       Y-dsp       Z-dsp       X-rot       Y-rot       Z-rot
     1  -3.665e-01  -8.573e-08  -7.767e-08   5.315e-09   8.214e-02   2.336e-09
     2  -1.214e-06  -1.816e-08   1.183e-06   2.610e-05  -6.070e-05   2.574e-05
     3  -1.213e-06   1.816e-08  -1.183e-06  -2.609e-05  -6.067e-05   2.573e-05
     4  -1.213e-06  -1.815e-08  -1.182e-06   2.609e-05  -6.067e-05  -2.573e-05
     5  -1.213e-06   1.815e-08   1.183e-06  -2.609e-05  -6.067e-05  -2.573e-05
  MODE     4:   f= 3.119372 Hz,  T= 0.320577 sec
		X- modal participation factor =   6.0069e-07 
		Y- modal participation factor =   4.8094e-07 
		Z- modal participation factor =   5.9990e-01 
  Joint   X-dsp       Y-dsp       Z-dsp       X-rot       Y-rot       Z-rot
     1  -4.624e-07  -6.864e-07   2.343e+00   6.626e-08  -5.687e-08  -1.547e-10
     2   3.332e-07   2.499e-07   3.212e-06   6.394e-05  -8.525e-05  -3.144e-10
     3  -3.324e-07   2.494e-07   3.205e-06   6.379e-05   8.506e-05  -2.403e-09
     4  -3.324e-07  -2.493e-07   3.205e-06  -6.380e-05   8.506e-05  -3.429e-10
     5   3.325e-07  -2.493e-07   3.205e-06  -6.379e-05  -8.506e-05   2.423e-09
  MODE     5:   f= 4.020972 Hz,  T= 0.248696 sec
		X- modal participation factor =   3.1571e-06 
		Y- modal participation factor =   2.7999e-01 
		Z- modal participation factor =   2.2771e-05 
  Joint   X-dsp       Y-dsp       Z-dsp       X-rot       Y-rot       Z-rot
     1  -1.472e-06   2.962e+00  -3.523e-06  -1.118e-01  -2.475e-07   4.786e-08
     2   6.454e-07   1.048e-06   5.786e-06   8.352e-05  -1.516e-04   1.738e-05
     3  -6.413e-07   1.041e-06   5.750e-06   8.304e-05   1.507e-04  -1.725e-05
     4   6.418e-07   1.041e-06  -5.750e-06   8.303e-05  -1.507e-04  -1.724e-05
     5  -6.414e-07   1.041e-06  -5.751e-06   8.305e-05   1.507e-04   1.724e-05
  MODE     6:   f= 4.706665 Hz,  T= 0.212465 sec
		X- modal participation factor =   4.1912e-01 
		Y- modal participation factor =   7.8173e-06 
		Z- modal participation factor =   5.6297e-05 
  Joint   X-dsp       Y-dsp       Z-dsp       X-rot       Y-rot       Z-rot
     1   2.847e+00  -3.120e-06  -6.201e-06   5.984e-07   7.390e-02  -2.481e-07
     2   4.416e-06   8.767e-07   7.529e-06   1.420e-04  -9.535e-05  -7.134e-05
     3   4.362e-06  -8.670e-07  -7.443e-06  -1.404e-04  -9.434e-05  -7.043e-05
     4   4.362e-06   8.674e-07  -7.441e-06   1.403e-04  -9.429e-05   7.043e-05
     5   4.362e-06  -8.662e-07   7.442e-06  -1.404e-04  -9.431e-05   7.046e-05
M A T R I X    I T E R A T I O N S: 2
There are 6 modes below 4.706665 Hz. ... All 6 modes were found.
