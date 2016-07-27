#!/usr/bin/python

# deltaE.py
# T.Trevethan
# 2016
#
# Script to calculate potential energy changes due to local strain from supplied 
# dipole tensors. Two 2D dipole tensors are supplied: gmin for changes to minima on the 
# PES and gbar for changes to saddle point energies. 
#
# The values given relate to the PES for the lattice vacancy in graphene, determined from
# GGA DFT models
#
# The atomistic strain tensors are read-in from the tensr.out file (created by gtensor.py)
# The script then outputs the modifications (in eV) to the minima for each 3 direction, and the 
# modifications to the saddle-point in each three direction, for each position. 

import sys
import math 
import argparse
import numpy as np

# define the mimima dipole tensor
gmin = np.array([
        [ -30.37 ,   0.0  ],
        [   0.0  ,   8.11 ]])

#define the saddle-point tensor
gbar = np.array([
        [  19.86 ,   0.0  ],
        [   0.0  , -51.26 ]])

print "gtensor: strain tensors to dE"

# Open input file
input  = open("tensr.out", 'r')
# Open output file
output = open("de.dat",'w')

#read input file and store tensors
tnsrs = []

line = input.readline()
#loop over all lines in input file
while line != "":
    data = line.split()
    tnsr = [float(data[0]),float(data[1]),float(data[2]),float(data[3]),float(data[4]),float(data[5]),float(data[6]),float(data[7]),float(data[8])]
    tnsrs.append(tnsr)
    line = input.readline()

print "Read in "+str(len(tnsrs))+" tensors"

for i in range(len(tnsrs)):
    emd1 = -tnsrs[i][0]*gmin[0,0] - tnsrs[i][1]*gmin[1,1] - tnsrs[i][2]*gmin[0,1]
    emd2 = -tnsrs[i][3]*gmin[0,0] - tnsrs[i][4]*gmin[1,1] - tnsrs[i][5]*gmin[0,1]
    emd3 = -tnsrs[i][6]*gmin[0,0] - tnsrs[i][7]*gmin[1,1] - tnsrs[i][8]*gmin[0,1]

    ebd1 = -tnsrs[i][0]*gbar[0,0] - tnsrs[i][1]*gbar[1,1] - tnsrs[i][2]*gbar[0,1]
    ebd2 = -tnsrs[i][3]*gbar[0,0] - tnsrs[i][4]*gbar[1,1] - tnsrs[i][5]*gbar[0,1]
    ebd3 = -tnsrs[i][6]*gbar[0,0] - tnsrs[i][7]*gbar[1,1] - tnsrs[i][8]*gbar[0,1]

    outline = str(emd1)+' '+str(emd2)+' '+str(emd3)+' '+str(ebd1)+' '+str(ebd2)+' '+str(ebd3)+"\n"
    output.write(outline)

print "Complete"
