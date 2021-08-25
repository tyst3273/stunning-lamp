# system modules
import numpy as np

# custom modules
import mod_lattice 
import mod_orbitals


# --------------------------------------------------------------------------------------------------
# set up the primitive cell of the crystal
# --------------------------------------------------------------------------------------------------

lattice = mod_lattice.lattice()

prim_scale = 3.929 # scaling factor for primitive vectors 
prim_vecs = [[1.000,0.000,0.000], # row vectors in angstroms
             [0.000,1.000,0.000],
             [0.000,0.000,2.496]]
lattice.set_primitive_vectors(prim_vecs,prim_scale) # set the primitive lattice vectors

# basis positions in crystal coords
basis_pos = [[0.000,0.000,0.000], # Cu
             [0.500,0.000,0.000], # O
             [0.000,0.500,0.000]] # O
basis_labels = ['Cu','O','O'] # optional metadata about the atoms
lattice.set_basis_positions(basis_pos,basis_labels)


# --------------------------------------------------------------------------------------------------
# set up the supercell that includes nearest neighbors
# --------------------------------------------------------------------------------------------------

sc_scale = 3.929 # scaling factor for supercell
sc_vecs = [[2.000,0.000,0.000], # row vectors angstroms
           [0.000,2.000,0.000],
           [0.000,0.000,2.496]]
lattice.set_supercell_vectors(sc_vecs,sc_scale) # set the supercell lattice vectors

# supercell positions in crystal coords
sc_pos = [[0.000,0.000,0.000], # Cu  
          [0.250,0.000,0.000], # O
          [0.000,0.250,0.000], # O
          [0.500,0.000,0.000], # Cu
          [0.750,0.000,0.000], # O
          [0.500,0.250,0.000], # O
          [0.500,0.500,0.000], # Cu
          [0.750,0.000,0.000], # O
          [0.500,0.750,0.000], # O
          [0.000,0.500,0.000], # Cu
          [0.250,0.500,0.000], # O
          [0.000,0.750,0.000]] # O
sc_labels = ['Cu','O','O','Cu','O','O','Cu','O','O','Cu','O','O'] # optional metadata about the atoms
lattice.set_supercell_positions(sc_pos,sc_labels)


# --------------------------------------------------------------------------------------------------
# find the neighbor lists in the supercell
# --------------------------------------------------------------------------------------------------

lattice.find_neighbors() # find neighbor lists using minimum image convention
lattice.write_neighbor_lists() # write neighbor lists to a file









