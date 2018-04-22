#!/usr/bin/python3

################################################################################################################
# FORMULAE FROM: https://www.supermagnete.de/eng/faq/How-do-you-calculate-the-magnetic-flux-density
################################################################################################################
'''
Formula for block magnet flux density
Formula for the B field on the symmetry axis of an axially magnetised block magnet (block or cube):


Br: Remanence field, independent of the magnet's geometry (see physical magnet data)
z: Distance from a pole face on the symmetry axis
L: Length of the block
W: Width of the block
D: Thickness (or height) of the block
The unit of length can be selected arbitrarily, as long as it is the same for all lengths.

################################################################################################################

Formula for cylinder magnet flux density
Formula for the B field on the symmetry axis of an axially magnetised cylinder magnet (disc or rod):

Br: Remanence field, independent of the magnet's geometry (see physical magnet data)
z: Distance from a pole face on the symmetrical axis
D: Thickness (or height) of the cylinder
R: Semi-diameter (radius) of the cylinder
The unit of length can be selected arbitrarily, as long as it is the same for all lengths.

'''

from math import sqrt, pi, atan
print("CYLINDER")
D = 0.003#*19
R = 0.01961504524593303 #0.0045
z = eval(input("enter distance in meters: "))

# https://www.nve.com/newsletter/nve_newsletter--7-12.htm
# Br:   Residual magnetism
Br = 1.17 # millitesla # N35 11700-12100 Gauss  1.17-1.21 Tesla 10.8-11.5 860-915 ≥12 ≥955  33-35 263-279 ≤80
# coercise magnet field for 1.17 T  is 931055 Am^1
B = (Br/2) * ( (D+z)/sqrt(R**2+(D+z)**2)  - (z/sqrt(R**2 + z**2) ))
print("CYLINDER magnet R", R,":\n\t\t\t", B, "T")

L = W = sqrt(pi * R**2)
B = Br/pi * ( atan(L * W/(2 * z * sqrt(4*z**2 + L**2 + W**2))) - atan(L * W/(2 * (D + z) * sqrt(4*(D + z)**2 + L**2 + W**2))) )
print("BLOCK magnet L",L,"W",W,":\n\t\t\t",B, "T")
