* Model of a single member for use in checking orientation commands
VERS     5
TYPE     5
VERT     3
UNIT     1 m   kN  T   C
 
NODE 1 -1 0.1 0    000000
NODE 2 -2 0.6 0.4  000000
NODE 3 -1 0.1 0    000000
NODE 4 -3 0.9 0.6  000000

MEMB 1 1 2 +Z 1 1 000000 000000
MEMB 2 3 4 +Z 1 1 000000 000000

MOFF 1 LO 0.1 0.1 0.1 0.1 0.1 0.1

PROP     1 PRIS LH80603010x1.2  Lipped Hat G350
       2.9180E-04   0.000       0.000       0.000      2.8270E-07  4.0540E-07

END

