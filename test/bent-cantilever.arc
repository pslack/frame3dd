* Bent cantilever
* Model created as test case for Microstran Parser in FRAME3DD
 
VERS     5
TYPE     5
VERT     3
UNIT     1 m   kN  T   C
 
NODE     1     0.0000     0.0000     0.0000 111111
NODE     2     0.0000     0.0000     1.0000 000000
NODE     3     0.0000     0.5000     1.5000 000000
 
MEMB     2     1     2    X      1     1    000000 000000
MEMB     3     3     2    X      1     1    000000 000000
 
PROP     1 LIBR asw       150UB14.0       Y default
      1.7800E-03  0.000      0.000     2.8100E-08 4.9500E-07 6.6600E-06
 
MATL     1 2.000E+08 2.500E-01 7.850E+00 1.170E-05
 
CASE     1 Lateral loading
NDLD     3      1.000      0.000     -1.000      0.000      0.000      0.000
 
END
