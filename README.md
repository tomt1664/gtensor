# gtensor

Script to calculate the 2D atomistic strain tensors for a supplied
sp2 bonded network. The numpy.linalg routines are used for the tensor 
calculation. 

The structure is supplied as an xyz format file: input.xyz

A bond table is constructed based on a maxmum bond length. This is set to 1.8 A by default, but can be over-ridden with the command line argument: 

       -c float : set the C-C bond length

Additionally, optional periodic boundary conditions can be set

       -x float : peridic boundary conditions in x-direction
       -y float : peridic boundary conditions in y-direction
       -z float : peridic boundary conditions in z-direction
 
The tensor is calculated relative to a perfect lattice structure. For an sp2 bonded system with three-coordinated atoms, this is dependent only on a single lattice parameter: in this case, the nearest neighbour bond distance. This is set to a default of 1.42 A (graphene, GGA DFT) but can be over-ridden with the command line argument:

       -b float : lattice bond distance
 
The three non-equivalent components of the 2D tensor (exx, eyy and exy (=eyx)) are then written to the file tensr.out for each of the three equivalent lattice directions in graphene. This gives 9 columns in total, with each row corresponding to the atomic positions. If an atom is not 3-coordinated, the entry is set to 99.0

> deltaE.py

Script to calculate potential energy changes due to local strain from supplied dipole tensors. Two 2D dipole tensors are supplied: gmin for changes to minima on the PES and gbar for changes to saddle point energies. 

The values given relate to the PES for the lattice vacancy in graphene, determined from GGA DFT models

The atomistic strain tensors are read-in from the tensr.out file (created by gtensor.py). The script then outputs the modifications (in eV) to the minima for each 3 direction, and the modifications to the saddle-point in each three direction, for each position. 
