* Set up frame
* test node forces for balance
 
VERS     5
TYPE     5
VERT     3
UNIT     1 m   kN  T   C
 
NODE     1     1.0000     0.5774     1.6334 000000
NODE     2     1.0000     1.7334     0.0000 111000
NODE     3     0.0000     0.0000     0.0000 111111
NODE     4     2.0000     0.0000     0.0000 111111
 
MEMB     2     2     3    Z      1     1    000000 000000
MEMB     3     3     4    Z      1     1    000000 000000
MEMB     4     4     2    Z      1     1    000000 000000
MEMB     5     1     3    Z      1     1    000000 000000
MEMB     6     1     4    Z      1     1    000000 000000
 
PROP     1 LIBR asw       60.3X2.3CHS     Y default
      4.1900E-04  0.000      0.000     3.5300E-07 1.7700E-07 1.7700E-07
 
MATL     1 2.000E+08 2.500E-01 7.850E+00 1.170E-05
 
CASE     1 Gravity
GRAV           0.000      0.000     -9.810
 
CASE     2 Nodal load
NDLD     1     10.000     40.000     20.000      0.000      0.000      0.000
 
CASE     3 Case 1 + Case2
COMB     1      1.000
COMB     2      1.000
 
END
 
$ design data
SMEM 2 CODE AS4100 DSEC CHS GR C350 MAXSL 180 DMAX 9999.0
     OFFS 0.000     TF L   BF L  LH S   KX 1.000   KY 1.000
     OFFS L         TF L   BF L  LH S
SMEM 3 CODE AS4100 DSEC CHS GR C350 MAXSL 180 DMAX 9999.0
     OFFS 0.000     TF L   BF L  LH S   KX 1.000   KY 1.000
     OFFS L         TF L   BF L  LH S
SMEM 4 CODE AS4100 DSEC CHS GR C350 MAXSL 180 DMAX 9999.0
     OFFS 0.000     TF L   BF L  LH S   KX 1.000   KY 1.000
     OFFS L         TF L   BF L  LH S
SMEM 5 CODE AS4100 DSEC CHS GR C350 MAXSL 180 DMAX 9999.0
     OFFS 0.000     TF L   BF L  LH S   KX 1.000   KY 1.000
     OFFS L         TF L   BF L  LH S
SMEM 6 CODE AS4100 DSEC CHS GR C350 MAXSL 180 DMAX 9999.0
     OFFS 0.000     TF L   BF L  LH S   KX 1.000   KY 1.000
     OFFS L         TF L   BF L  LH S
END
