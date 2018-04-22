#!/bin/bash
#awk '$3==-0.4 { print "" } 1' bm_box.dat.csv > bm_box.dat.csv.4d
#awk '(sqrt((sqrt($3*$3+$4*$4+$5*$5) - 0.3)^2) < 0.02) || !NF' bm_box.dat.csv.4d > bm_box.dat.csv.4d.sphere
#gnuplot -e "set xrange[-0.4:0.4];set yrange[-0.4:0.4];set zrange[-0.4:0.4];set parametric;set hidden3d offset 0;set xlabel 'x';set ylabel 'y';set zlabel 'z';set view equal xyz;set term wxt;set palette defined (-1 'blue',0 'white', 1 'red');set palette maxcolors 25 ; set pm3d interpolate 0,0;splot  'bm_box.dat.csv.4d.sphere' using 3:4:5:(column(9)>0.00000000003?(column(9)):1/0) pt 7 ps .05  palette" -p
gnuplot -e "set xrange[-0.4:0.4];set yrange[-0.4:0.4];set zrange[-0.4:0.4];set parametric;set hidden3d offset 0;set xlabel 'x';set ylabel 'y';set zlabel 'z';set view equal xyz;set term wxt;set palette defined (-1 'blue',0 'white', 1 'red');set palette maxcolors 25 ; set pm3d interpolate 0,0;splot  'bm_box.dat.csv.4d.sphere' using 3:4:5:(column(9)>0.00000000009?(column(9)):1/0) pt 7 ps .05  palette" -p
