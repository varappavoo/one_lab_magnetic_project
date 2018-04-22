#!/bin/bash
#awk '$3==-1 { print "" } 1' b_box.dat.csv > b_box.dat.csv.4d
#awk '(sqrt((sqrt($3*$3+$4*$4+$5*$5) - 0.35)^2) < 0.02) || !NF' b_box.dat.csv.4d > b_box.dat.csv.4d.sphere
#awk '{print $3,$4,$5,sqrt($9*$9+$10*$10+$11*$11)}' b_box.dat.csv.4d.sphere > b_box.dat.csv.4d.sphere.mag
gnuplot -e ";set cbtics format '%g';;set xrange[-0.4:0.4];set yrange[-0.4:0.4];set zrange[-0.4:0.4];set parametric;set hidden3d offset 0;set xlabel 'x';set ylabel 'y';set zlabel 'z';set view equal xyz;set term wxt; set palette defined (-1 'white',0 'red', 1 'black');;set palette maxcolors 25 ; set pm3d interpolate 0,0;splot  'b_box.dat.csv.4d.sphere.mag' every 25 using 1:2:3:(column(4)>0.0?(column(4)):1/0) pt 7 ps .1  palette" -p

# 2D VIEW
gnuplot -e "set cbrange [0.00001:];set cbtics format '%g';set xrange[-0.4:0.4];set yrange[-0.4:0.4];set zrange[-0.4:0.4];set hidden3d offset 0;set xlabel 'x';set ylabel 'y';set zlabel 'z';set view equal xyz; set palette defined (-1 'white',0 'red', 1 'black');;set palette maxcolors 25 ; set pm3d interpolate 0,0;plot  'b_box.dat.csv.4d.sphere.mag' every 1 using 1:2:(column(4)>0.00001?(column(4)):1/0) pt 7 ps .5  palette" -p
gnuplot -e "set cbrange [0.00001:];set cbtics format '%g';set xrange[-0.4:0.4];set yrange[-0.4:0.4];set zrange[-0.4:0.4];set hidden3d offset 0;set xlabel 'x';set ylabel 'z';set zlabel 'z';set view equal xyz; set palette defined (-1 'white',0 'red', 1 'black');;set palette maxcolors 25 ; set pm3d interpolate 0,0;plot  'b_box.dat.csv.4d.sphere.mag' every 1 using 1:3:(column(4)>0.00001?(column(4)):1/0) pt 7 ps .5  palette" -p
gnuplot -e "set cbrange [0.00001:];set cbtics format '%g';set xrange[-0.4:0.4];set yrange[-0.4:0.4];set zrange[-0.4:0.4];set hidden3d offset 0;set xlabel 'y';set ylabel 'z'; set palette defined (-1 'white',0 'red', 1 'black');;set palette maxcolors 25 ; set pm3d interpolate 0,0;plot  'b_box.dat.csv.4d.sphere.mag' every 1 using 2:3:(column(4)>0.00001?(column(4)):1/0) pt 7 ps .5  palette" -p



#gnuplot -e "set xrange[-0.4:0.4];set yrange[-0.4:0.4];set zrange[-0.4:0.4];set parametric;set hidden3d offset 0;set xlabel 'x';set ylabel 'y';set zlabel 'z';set view equal xyz;set term wxt; set palette defined (-1 'yellow',0 'red', 1 'black');;set palette maxcolors 25 ; set pm3d interpolate 0,0;splot  'bm_box.dat.csv.4d.sphere.mag' using 3:4:5:(column(9)>0.00000000003?(column(9)):1/0) pt 7 ps .1  palette" -p
#vector plot
#awk '$3==-1 { print "" } 1' b_box.dat.csv > b_box.dat.csv.4d
#gnuplot -e "set xrange[-0.4:0.4];set yrange[-0.4:0.4];set zrange[-0.4:0.4];set parametric;set hidden3d offset 0;set xlabel 'x';set ylabel 'y';set zlabel 'z';set view equal xyz; set palette defined (-1 'yellow',0 'red', 1 'black');set palette maxcolors 25 ; set pm3d interpolate 0,0;set term wxt;splot  'b_box.dat.csv.4d' using 3:4:5:9:10:(sqrt((column(9)**2)+(column(10)**2)+(column(11)**2))<0.000001?(column(11)):1/0) with vectors filled head lw 3 lc palette" -p
#gnuplot -e "sample=20;scale=1000;set xrange[-0.4:0.4];set yrange[-0.4:0.4];set zrange[-0.4:0.4];set parametric;set hidden3d offset 0;set xlabel 'x';set ylabel 'y';set zlabel 'z';set view equal xyz; set palette defined (-1 'yellow',0 'red', 1 'black');set palette maxcolors 25 ; set pm3d interpolate 0,0;set term wxt;splot  'b_box.dat.csv.4d' every sample using 3:4:5:(column(9)*scale):(column(10)*scale):(sqrt((column(9)**2)+(column(10)**2)+(column(11)**2))>0.000001?(column(11)*scale):1/0) with vectors filled head lw 1 lc palette" -p
#gnuplot -e "set cbtics format '%g';set term wxt;set grid xtics, ytics, ztics; set size 1,1;set xrange[-0.4:0.4];set yrange[-0.4:0.4];set zrange[-0.4:0.4];set parametric;set hidden3d offset 10;set xlabel 'x';set ylabel 'y';set zlabel 'z';set view equal xyz; set palette defined (-1 'yellow',0 'red', 1 'black');set palette maxcolors 25 ; set pm3d interpolate 0,0;splot  'bm_box.dat.csv.4d.sphere' using 3:4:5:(column(9)>0.00000000003?(column(9)):1/0) pt 7 ps .3 palette" -p


