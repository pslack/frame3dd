# FRAME3DD ANALYSIS RESULTS  http://frame3dd.sf.net/
# frame analysis via Matlab interface 
# Mon Apr  8 05:02:00 1974
# M E S H   A N N O T A T I O N   F I L E 
set title "frame analysis via Matlab interface \nanalysis file: IOdata.FMM   deflection exaggeration: 10.0   load case 1 of 1 "
set autoscale
set noborder
set pointsize 1.0
set xtics; set ytics; set ztics; 
set nozeroaxis
set nokey
set nolabel
# NODE NUMBER LABELS
set label ' 1' at   0.0000e+00,   0.0000e+00,   0.0000e+00
set label ' 2' at   1.2000e+01,   0.0000e+00,   0.0000e+00
set label ' 3' at   2.4000e+01,   0.0000e+00,   0.0000e+00
set label ' 4' at   3.6000e+01,   0.0000e+00,   0.0000e+00
set label ' 5' at   4.8000e+01,   0.0000e+00,   0.0000e+00
set label ' 6' at   6.0000e+01,   0.0000e+00,   0.0000e+00
set label ' 7' at   7.2000e+01,   0.0000e+00,   0.0000e+00
set label ' 8' at   1.2000e+01,   1.2000e+01,   0.0000e+00
set label ' 9' at   2.4000e+01,   1.2000e+01,   0.0000e+00
set label ' 10' at   3.6000e+01,   1.2000e+01,   0.0000e+00
set label ' 11' at   4.8000e+01,   1.2000e+01,   0.0000e+00
set label ' 12' at   6.0000e+01,   1.2000e+01,   0.0000e+00
# MEMBER NUMBER LABELS
set label ' 1' at   6.0000e+00,   0.0000e+00,   0.0000e+00
set label ' 2' at   1.8000e+01,   0.0000e+00,   0.0000e+00
set label ' 3' at   3.0000e+01,   0.0000e+00,   0.0000e+00
set label ' 4' at   4.2000e+01,   0.0000e+00,   0.0000e+00
set label ' 5' at   5.4000e+01,   0.0000e+00,   0.0000e+00
set label ' 6' at   6.6000e+01,   0.0000e+00,   0.0000e+00
set label ' 7' at   6.0000e+00,   6.0000e+00,   0.0000e+00
set label ' 8' at   1.2000e+01,   6.0000e+00,   0.0000e+00
set label ' 9' at   1.8000e+01,   6.0000e+00,   0.0000e+00
set label ' 10' at   2.4000e+01,   6.0000e+00,   0.0000e+00
set label ' 11' at   3.0000e+01,   6.0000e+00,   0.0000e+00
set label ' 12' at   3.6000e+01,   6.0000e+00,   0.0000e+00
set label ' 13' at   4.2000e+01,   6.0000e+00,   0.0000e+00
set label ' 14' at   4.8000e+01,   6.0000e+00,   0.0000e+00
set label ' 15' at   5.4000e+01,   6.0000e+00,   0.0000e+00
set label ' 16' at   6.0000e+01,   6.0000e+00,   0.0000e+00
set label ' 17' at   6.6000e+01,   6.0000e+00,   0.0000e+00
set label ' 18' at   1.8000e+01,   1.2000e+01,   0.0000e+00
set label ' 19' at   3.0000e+01,   1.2000e+01,   0.0000e+00
set label ' 20' at   4.2000e+01,   1.2000e+01,   0.0000e+00
set label ' 21' at   5.4000e+01,   1.2000e+01,   0.0000e+00
plot '/tmp/IOdata-msh' u 2:3 t 'undeformed mesh' w lp lw 1 lt 5 pt 6, '/tmp/IOdata-mshf.001' u 1:2 t 'load case 1 of 1' w l lw 2 lt 3
# set parametric
# set view 60, 70, 1 
# set nokey
# set xlabel 'x'
# set ylabel 'y'
# set zlabel 'z'
# splot '/tmp/IOdata-msh' u 2:3:4 t 'load case 1 of 1' w lp  lw 1 lt 5 pt 6, '/tmp/IOdata-mshf.001' u 1:2:3 t 'load case 1 of 1' w l lw 2 lt 3
