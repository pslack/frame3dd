* Lateral column
* Simple test case for demonstrating Microstran parser in Frame3d
 
VERS     5
TYPE     5
VERT     3
UNIT     1 m   kN  T   C
 
NODE     1     0.0000     0.0000     0.0000 111111
NODE     2    -1.0000     0.0000     0.0000 000000
 
MEMB     1     1     2   -Z      1     1    000000 000000
 
PROP     1 LIBR asw       65X35X2.0RHS    Y comment
      3.7400E-04  0.000      0.000     1.8400E-07 7.7800E-08 2.0400E-07
 
MATL     1 2.000E+08 2.500E-01 7.850E+00 1.170E-05
 
CASE     1 Horizontal load 1kN
NDLD     2     -3.000     -3.000      5.000      0.000      0.000      0.000
 
END
