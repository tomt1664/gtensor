#!/usr/bin/python

# gtensor.py
# T.Trevethan
# 2016
#
# Python code to calculate the 2D atomistic strain tensors for a supplied
# sp2 bonded network. The numpy.linalg routines are used for the tensor 
# calculation. 
#
# The structure is supplied as an xyz format file: input.xyz
# 
# A bond table is constructed based on a maxmum bond length. This is set to
# 1.8 A by default, but can be over-ridden with the command line argument: 
#
#       -c float : set the C-C bond length
#
# Additionally, optional periodic boundary conditions can be set
#
#       -x float : peridic boundary conditions in x-direction
#       -y float : peridic boundary conditions in y-direction
#       -z float : peridic boundary conditions in z-direction
# 
# The tensor is calculated relative to a perfect lattice structure. For an sp2 
# bonded with three-coordinated atoms, this is dependent only on a single lattice 
# parameter: in this case, the nearest neighbour bond distance. This is set to a default
# of 1.42 A (graphene, GGA DFT) but can be over-ridden with the command line argument:
#
#       -b float : lattice bond distance
# 
# The three non-equivalent components of the 2D tensor (exx, eyy and exy (=eyx)) are then 
# written to the file tensr.out for each of the three equivalent lattice directions in graphene. 
# This gives 9 columns in total, with each row corresponding to the atomic positions. 
#
# If an atom is not 3-coordinated, the entry is set to 99.0

import sys
import math 
import argparse
import numpy as np

#function to rotate a point about the origin in the x-y plane
def rotate(xcoord,ycoord,ang):
    xrot = xcoord*math.cos(ang) - ycoord*math.sin(ang)
    yrot = xcoord*math.sin(ang) + ycoord*math.cos(ang)
    rot = [xrot,yrot]
    return rot

cutoff = 1.7
latbnd = 1.42
pi = 3.141592654
xper = 0
yper = 0
zper = 0

# Parse command line input
parser = argparse.ArgumentParser()
parser.add_argument('-b', type = float, help = 'lattice bond length')
parser.add_argument('-c', type = float, help = 'bond cutoff')
parser.add_argument('-x', type = float, help = 'x-direction PBC')
parser.add_argument('-y', type = float, help = 'y-direction PBC')
parser.add_argument('-z', type = float, help = 'z-direction PBC')
variables = parser.parse_args()
if variables.b:
    latbnd = variables.b
if variables.c:
    cutoff = variables.c
if variables.x:
    xper = variables.x
if variables.y:
    yper = variables.y
if variables.z:
    zper = variables.z

print "gtensor: 2D atomistic strain tensors"

per = (xper,yper,zper)

# Open input file
input  = open("input.xyz", 'r')
# Open output file
output = open("tensr.out",'w')

#read input file and store carbon atom coordinates in array
coords = []
label = []

numline = input.readline()
nat = int(numline)
blank = input.readline()

# read coordinates
nc = 0
for i in range(nat):
    line = input.readline()
    fdata = line.split()
    data = (float(fdata[1]),float(fdata[2]),float(fdata[3]))
    coords.append(data)
    label.append(fdata[0])

print "Read in "+str(nat)+" atoms"

# define the lattice coordination structures
ulcs = []
tdat = (0.0,latbnd)
ulcs.append(tdat)
tdat = (latbnd*math.sin(pi/3),-latbnd*0.5)
ulcs.append(tdat)
tdat = (-latbnd*math.sin(pi/3),-latbnd*0.5)
ulcs.append(tdat)

#create bond table
bond = []

for i in range(nat):
    ib = 0
    bond.append([])
    for j in range(nat):
        if i != j:
            p1 = [coords[i][0],coords[i][1],coords[i][2]]
            p2 = [coords[j][0],coords[j][1],coords[j][2]]
# apply boundary conditions
            for n in range(2):
                if per[n] > 0:
                    if p1[n] < per[n]/4 and p2[n] > 3*per[n]/4:
                        p2[n]-= per[n]
                    if p2[n] < per[n]/4 and p1[n] > 3*per[n]/4:
                        p1[n] -= per[n]
            xdst = p1[0] - p2[0]
            ydst = p1[1] - p2[1]
            zdst = p1[2] - p2[2]
            dst = math.sqrt(xdst*xdst+ydst*ydst+zdst*zdst)
            if dst < cutoff:
                if ib > 2:
                    print "Error: atom "+str(i)+" has more than 3 bonds"
                    sys.exit(1)
                bond[i].append(j)
                ib += 1


# calculate strain tensors
for i in range(nat):
    if len(bond[i]) == 3:
# define three bond vectors
        p1 = [coords[i][0],coords[i][1]]
        p2 = [coords[bond[i][0]][0],coords[bond[i][0]][1]]
        p3 = [coords[bond[i][1]][0],coords[bond[i][1]][1]]
        p4 = [coords[bond[i][2]][0],coords[bond[i][2]][1]]
# apply boundary conditions
        for n in range(2):
            if per[n] > 0:
                if p1[n] < per[n]/4 and p2[n] > 3*per[n]/4:
                    p2[n]-= per[n]
                if p2[n] < per[n]/4 and p1[n] > 3*per[n]/4:
                    p1[n] -= per[n]
                if p1[n] < per[n]/4 and p3[n] > 3*per[n]/4:
                    p3[n]-= per[n]
                if p3[n] < per[n]/4 and p1[n] > 3*per[n]/4:
                    p1[n] -= per[n]
                if p1[n] < per[n]/4 and p4[n] > 3*per[n]/4:
                    p4[n]-= per[n]
                if p4[n] < per[n]/4 and p1[n] > 3*per[n]/4:
                    p1[n] -= per[n]
        vec = []
        vec1 = [p2[0] - p1[0], p2[1] - p1[1]]
        vec2 = [p3[0] - p1[0], p3[1] - p1[1]]
        vec3 = [p4[0] - p1[0], p4[1] - p1[1]]
        vec.append(vec1)
        vec.append(vec2)
        vec.append(vec3)
# find the principle cooridinator and orientation
        if vec1[1] > 0.0 and vec2[1] < 0.0 and vec3[1] < 0.0:
            ip1 = 0
            ud = 1
            if vec2[0] > 0.0 and vec3[0] < 0.0:
                ip2 = 1
                ip3 = 2
            elif vec3[0] > 0.0 and vec2[0] < 0.0:
                ip2 = 2
                ip3 = 1
        elif vec2[1] > 0.0 and vec1[1] < 0.0 and vec3[1] < 0.0:
            ip1 = 1
            ud = 1
            if vec1[0] > 0.0 and vec3[0] < 0.0:
                ip2 = 0
                ip3 = 2
            elif vec3[0] > 0.0 and vec1[0] < 0.0:
                ip2 = 2
                ip3 = 0
        elif vec3[1] > 0.0 and vec1[1] < 0.0 and vec2[1] < 0.0:
            ip1 = 2
            ud = 1
            if vec1[0] > 0.0 and vec2[0] < 0.0:
                ip2 = 0
                ip3 = 1
            elif vec2[0] > 0.0 and vec1[0] < 0.0:
                ip2 = 1
                ip3 = 0
        elif vec1[1] < 0.0 and vec2[1] > 0.0 and vec3[1] > 0.0:
            ip1 = 0
            ud = -1
            if vec2[0] > 0.0 and vec3[0] < 0.0:
                ip2 = 1
                ip3 = 2
            elif vec3[0] > 0.0 and vec2[0] < 0.0:
                ip2 = 2
                ip3 = 1
        elif vec2[1] < 0.0 and vec1[1] > 0.0 and vec3[1] > 0.0:
            ip1 = 1
            ud = -1
            if vec1[0] > 0.0 and vec3[0] < 0.0:
                ip2 = 0
                ip3 = 2
            elif vec3[0] > 0.0 and vec1[0] < 0.0:
                ip2 = 2
                ip3 = 0
        elif vec3[1] < 0.0 and vec1[1] > 0.0 and vec2[1] > 0.0:
            ip1 = 2
            ud = -1
            if vec1[0] > 0.0 and vec2[0] < 0.0:
                ip2 = 0
                ip3 = 1
            elif vec2[0] > 0.0 and vec1[0] < 0.0:
                ip2 = 1
                ip3 = 0
# set the lattice direction
        vec1[1] = vec1[1]*ud
        vec2[1] = vec2[1]*ud
        vec3[1] = vec3[1]*ud
# reduce coordinates to two lattice vectors
        latvec1 = [vec[ip1][0] - vec[ip3][0],vec[ip1][1] - vec[ip3][1]]
        latvec2 = [vec[ip2][0] - vec[ip3][0],vec[ip2][1] - vec[ip3][1]]
# solve system
        amat = np.array([
                [ulcs[0][0] - ulcs[2][0], 0, ulcs[1][0] - ulcs[2][0], 0],
                [ulcs[0][1] - ulcs[2][1], 0, ulcs[1][1] - ulcs[2][1], 0],
                [0, ulcs[0][0] - ulcs[2][0], 0, ulcs[1][0] - ulcs[2][0]],
                [0, ulcs[0][1] - ulcs[2][1], 0, ulcs[1][1] - ulcs[2][1]]])
        
        bvec = np.array([latvec1[0],latvec1[1],latvec2[0],latvec2[1]])
        cmat = np.transpose(amat)
        xvac = np.linalg.solve(cmat,bvec)
# save the tensor components: exx, eyy, exy
        tnsr1 = [xvac[0],xvac[3],xvac[1]*xvac[3] + xvac[2]*xvac[0] + xvac[1]*xvac[2]]

# rotate 60 degrees
# rotate lattice
        ulcsr = []
        ulcr1 = rotate(ulcs[0][0],ulcs[0][1],pi/3)
        ulcr2 = rotate(ulcs[1][0],ulcs[1][1],pi/3)
        ulcr3 = rotate(ulcs[2][0],ulcs[2][1],pi/3)
        ulcsr.append(ulcr1)
        ulcsr.append(ulcr2)
        ulcsr.append(ulcr3)
# rotate configuration
        vecr = []
        vec1r = rotate(vec1[0],vec1[1],pi/3)
        vec2r = rotate(vec2[0],vec2[1],pi/3)
        vec3r = rotate(vec3[0],vec3[1],pi/3)
        vecr.append(vec1r)
        vecr.append(vec2r)
        vecr.append(vec3r)
# reduce coordinates to two lattice vectors
        latvec1 = [vecr[ip1][0] - vecr[ip3][0],vecr[ip1][1] - vecr[ip3][1]]
        latvec2 = [vecr[ip2][0] - vecr[ip3][0],vecr[ip2][1] - vecr[ip3][1]]
# solve system
        amat = np.array([
                [ulcsr[0][0] - ulcsr[2][0], 0, ulcsr[1][0] - ulcsr[2][0], 0],
                [ulcsr[0][1] - ulcsr[2][1], 0, ulcsr[1][1] - ulcsr[2][1], 0],
                [0, ulcsr[0][0] - ulcsr[2][0], 0, ulcsr[1][0] - ulcsr[2][0]],
                [0, ulcsr[0][1] - ulcsr[2][1], 0, ulcsr[1][1] - ulcsr[2][1]]])
        
        bvec = np.array([latvec1[0],latvec1[1],latvec2[0],latvec2[1]])
        cmat = np.transpose(amat)
        xvac = np.linalg.solve(cmat,bvec)
# save the tensor components: exx, eyy, exy
        tnsr2 = [xvac[0],xvac[3],xvac[1]*xvac[3] + xvac[2]*xvac[0] + xvac[1]*xvac[2]]

# rotate 120 degrees
# rotate lattice
        ulcsr = []
        ulcr1 = rotate(ulcs[0][0],ulcs[0][1],2*pi/3)
        ulcr2 = rotate(ulcs[1][0],ulcs[1][1],2*pi/3)
        ulcr3 = rotate(ulcs[2][0],ulcs[2][1],2*pi/3)
        ulcsr.append(ulcr1)
        ulcsr.append(ulcr2)
        ulcsr.append(ulcr3)
# rotate configuration
        vecr = []
        vec1r = rotate(vec1[0],vec1[1],2*pi/3)
        vec2r = rotate(vec2[0],vec2[1],2*pi/3)
        vec3r = rotate(vec3[0],vec3[1],2*pi/3)
        vecr.append(vec1r)
        vecr.append(vec2r)
        vecr.append(vec3r)
# reduce coordinates to two lattice vectors
        latvec1 = [vecr[ip1][0] - vecr[ip3][0],vecr[ip1][1] - vecr[ip3][1]]
        latvec2 = [vecr[ip2][0] - vecr[ip3][0],vecr[ip2][1] - vecr[ip3][1]]
# solve system
        amat = np.array([
                [ulcsr[0][0] - ulcsr[2][0], 0, ulcsr[1][0] - ulcsr[2][0], 0],
                [ulcsr[0][1] - ulcsr[2][1], 0, ulcsr[1][1] - ulcsr[2][1], 0],
                [0, ulcsr[0][0] - ulcsr[2][0], 0, ulcsr[1][0] - ulcsr[2][0]],
                [0, ulcsr[0][1] - ulcsr[2][1], 0, ulcsr[1][1] - ulcsr[2][1]]])
        
        bvec = np.array([latvec1[0],latvec1[1],latvec2[0],latvec2[1]])
        cmat = np.transpose(amat)
        xvac = np.linalg.solve(cmat,bvec)
# save the tensor components: exx, eyy, exy
        tnsr3 = [xvac[0],xvac[3],xvac[1]*xvac[3] + xvac[2]*xvac[0] + xvac[1]*xvac[2]]
    else:
        tnsr1 = [99.0,99.0,99.0]
        tnsr2 = [99.0,99.0,99.0]
        tnsr3 = [99.0,99.0,99.0]

# write to file
    outline = ' '.join(map(str,tnsr1))+" "+' '.join(map(str,tnsr2))+" "+' '.join(map(str,tnsr3))+"\n"
    output.write(outline)

print "Complete"
