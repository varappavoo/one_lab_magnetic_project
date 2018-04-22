#!/bin/bash
awk '$3==-1 { print "" } 1' bm_box.dat.csv > bm_box.dat.csv.4d
awk '(sqrt((sqrt($3*$3+$4*$4+$5*$5) - 0.4)^2) < 0.02) || !NF' bm_box.dat.csv.4d > bm_box.dat.csv.4d.sphere
