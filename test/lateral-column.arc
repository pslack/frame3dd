* Lateral column
* Simple test case for demonstrating Microstran parser in Frame3dd
 
VERS     5
TYPE     5
VERT     3
UNIT     1 m   kN  T   C
 
NODE     1     0.0000     0.0000     0.0000 111111
NODE     2     0.0000     0.0000     1.0000 000000
 
MEMB     1     1     2    X      1     1    000000 000000
 
PROP     1 LIBR asw       50X20X1.6RHS    Y comment
      2.0700E-04  0.000      0.000     3.8900E-08 1.4200E-08 6.0800E-08
 
MATL     1 2.000E+08 2.500E-01 7.850E+00 1.170E-05
 
CASE     1 Horizontal load 1kN
NDLD     2      1.000      0.500      0.000      0.000      0.000      0.000
 
END
