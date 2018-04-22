#!/usr/bin/python3
# https://stackoverflow.com/questions/35402609/point-on-circle-base-on-given-angle


import sys
from math import cos, sin, pi

def point_on_circle(center, angle, radius):
    '''
        Finding the x,y coordinates on circle, based on given angle
    '''
    #center of circle, angle in degree and radius of circle
    x = center[0] + (radius * cos(angle))
    y = center[1] + (radius * sin(angle))
    return x,y

center = [float(sys.argv[1]), float(sys.argv[2])]
angle = (float(sys.argv[3])/180) * pi # degree to radian
radius = float(sys.argv[4])

print(point_on_circle(center, angle, radius))