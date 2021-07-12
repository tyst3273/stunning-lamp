import mod_lattice
import mod_forceconstants
import mod_dynmat 

# ---------------------------------------------------------------------------------------------------

"""
set up the lattice. need primitive cell vectors positions of atoms in primitive cell,
supercell vectors, positions in supercell, and set data for each atom type

"""

lattice = mod_lattice.lattice()

# set data for each atom type
atom_masses = [51.9961] # amu
lattice.set_atom_data(atom_masses)

# set primitive cell lattice vectors
prim_vecs = [[1,0,0], # angstrom
             [0,1,0],
             [0,0,1]]
prim_scale = 2.8525 # prim_vecs = prim_vecs*prim_scale
lattice.set_prim_vecs(prim_vecs,prim_scale)

# set basis positions in primitive cell
basis_pos = [[0.0,0.0,0.0], # crystal coords (prim_vecs)
             [0.5,0.5,0.5]]
basis_types = [1,1] # 1 for each atom
lattice.set_basis_pos(basis_pos,basis_types)

# set supercell lattice vectors
sc_vecs = [[1,0,0], # angstrom
           [0,1,0],
           [0,0,1]]
sc_scale = 5.7049
lattice.set_sc_vecs(sc_vecs,sc_scale)

# set positions of atoms in supercell
sc_pos = [[0.00,0.00,0.00], # crystal coords (sc_vecs)
          [0.50,0.00,0.00],
          [0.00,0.50,0.00],
          [0.50,0.50,0.00],
          [0.00,0.00,0.50],
          [0.50,0.00,0.50],
          [0.00,0.50,0.50],
          [0.50,0.50,0.50],
          [0.25,0.25,0.25],
          [0.75,0.25,0.25],
          [0.25,0.75,0.25],
          [0.75,0.75,0.25],
          [0.25,0.25,0.75],
          [0.75,0.25,0.75],
          [0.25,0.75,0.75],
          [0.75,0.75,0.75]]
sc_types = [1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1]
lattice.set_sc_pos(sc_pos,sc_types)

# -------------------------------------------------------------------------------------------------

"""
set force constants. there should be 1 between each primitive cell atom and supercell atom.

fc_ij = [[ phi(i,j)_x,x , phi(i,j)_x,y , phi(i,j)_x,z ],
         [ phi(i,j)_y,x , phi(i,j)_y,y , phi(i,j)_y,z ],
         [ phi(i,j)_z,x , phi(i,j)_z,y , phi(i,j)_z,z ]]

fc units = eV/(Angst*Angstr)

"""

fc_mats = mod_forceconstants.force_constant_matrices(lattice)


# set fc between atoms
basis_ind = 0 # which basis atom 
sc_ind_i = 0 # ind of basis atom in sc
sc_ind_j = 0 # ind of other atom coupled to basis atom
fc = [[ 12.065165, -0.000000, -0.000000],
      [ -0.000000, 12.065165, -0.000000],
      [ -0.000000, -0.000000, 12.065165]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 0
sc_ind_i = 0
sc_ind_j = 1
fc = [[ -5.806100, -0.000000, -0.000000],
      [ -0.000000,  0.402342, -0.000000],
      [ -0.000000, -0.000000,  0.402342]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 0
sc_ind_i = 0
sc_ind_j = 2
fc = [[  0.402342, -0.000000, -0.000000],
      [ -0.000000, -5.806100, -0.000000],
      [ -0.000000, -0.000000,  0.402342]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 0
sc_ind_i = 0
sc_ind_j = 3
fc = [[ -0.676651, -0.000000, -0.000000],
      [ -0.000000, -0.676651, -0.000000],
      [ -0.000000, -0.000000,  0.013278]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 0
sc_ind_i = 0
sc_ind_j = 4
fc = [[  0.402342, -0.000000, -0.000000],
      [ -0.000000,  0.402342, -0.000000],
      [ -0.000000, -0.000000, -5.806100]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 0
sc_ind_i = 0
sc_ind_j = 5
fc = [[ -0.676651, -0.000000, -0.000000],
      [ -0.000000,  0.013278, -0.000000],
      [ -0.000000, -0.000000, -0.676651]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 0
sc_ind_i = 0
sc_ind_j = 6
fc = [[  0.013278, -0.000000, -0.000000],
      [ -0.000000, -0.676651, -0.000000],
      [ -0.000000, -0.000000, -0.676651]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 0
sc_ind_i = 0
sc_ind_j = 7
fc = [[ -0.407391, -0.000000, -0.000000],
      [ -0.000000, -0.407391, -0.000000],
      [ -0.000000, -0.000000, -0.407391]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 0
sc_ind_i = 0
sc_ind_j = 8
fc = [[ -0.664542, -0.332266, -0.332266],
      [ -0.332266, -0.664542, -0.332266],
      [ -0.332266, -0.332266, -0.664542]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 0
sc_ind_i = 0
sc_ind_j = 9
fc = [[ -0.664542,  0.332266,  0.332266],
      [  0.332266, -0.664542, -0.332266],
      [  0.332266, -0.332266, -0.664542]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 0
sc_ind_i = 0
sc_ind_j = 10
fc = [[ -0.664542,  0.332266, -0.332266],
      [  0.332266, -0.664542,  0.332266],
      [ -0.332266,  0.332266, -0.664542]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 0
sc_ind_i = 0
sc_ind_j = 11
fc = [[ -0.664542, -0.332266,  0.332266],
      [ -0.332266, -0.664542,  0.332266],
      [  0.332266,  0.332266, -0.664542]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 0
sc_ind_i = 0
sc_ind_j = 12
fc = [[ -0.664542, -0.332266,  0.332266],
      [ -0.332266, -0.664542,  0.332266],
      [  0.332266,  0.332266, -0.664542]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 0
sc_ind_i = 0
sc_ind_j = 13
fc = [[ -0.664542,  0.332266, -0.332266],
      [  0.332266, -0.664542,  0.332266],
      [ -0.332266,  0.332266, -0.664542]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 0
sc_ind_i = 0
sc_ind_j = 14
fc = [[ -0.664542,  0.332266,  0.332266],
      [  0.332266, -0.664542, -0.332266],
      [  0.332266, -0.332266, -0.664542]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 0
sc_ind_i = 0
sc_ind_j = 15
fc = [[ -0.664542, -0.332266, -0.332266],
      [ -0.332266, -0.664542, -0.332266],
      [ -0.332266, -0.332266, -0.664542]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 1
sc_ind_i = 8
sc_ind_j = 0
fc = [[ -0.664542, -0.332266, -0.332266],
      [ -0.332266, -0.664542, -0.332266],
      [ -0.332266, -0.332266, -0.664542]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 1
sc_ind_i = 8
sc_ind_j = 1
fc = [[ -0.664542,  0.332266,  0.332266],
      [  0.332266, -0.664542, -0.332266],
      [  0.332266, -0.332266, -0.664542]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 1
sc_ind_i = 8
sc_ind_j = 2
fc = [[ -0.664542,  0.332266, -0.332266],
      [  0.332266, -0.664542,  0.332266],
      [ -0.332266,  0.332266, -0.664542]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 1
sc_ind_i = 8
sc_ind_j = 3
fc = [[ -0.664542, -0.332266,  0.332266],
      [ -0.332266, -0.664542,  0.332266],
      [  0.332266,  0.332266, -0.664542]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 1
sc_ind_i = 8
sc_ind_j = 4
fc = [[ -0.664542, -0.332266,  0.332266],
      [ -0.332266, -0.664542,  0.332266],
      [  0.332266,  0.332266, -0.664542]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 1
sc_ind_i = 8
sc_ind_j = 5
fc = [[ -0.664542,  0.332266, -0.332266],
      [  0.332266, -0.664542,  0.332266],
      [ -0.332266,  0.332266, -0.664542]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 1
sc_ind_i = 8
sc_ind_j = 6
fc = [[ -0.664542,  0.332266,  0.332266],
      [  0.332266, -0.664542, -0.332266],
      [  0.332266, -0.332266, -0.664542]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 1
sc_ind_i = 8
sc_ind_j = 7
fc = [[ -0.664542, -0.332266, -0.332266],
      [ -0.332266, -0.664542, -0.332266],
      [ -0.332266, -0.332266, -0.664542]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 1
sc_ind_i = 8
sc_ind_j = 8
fc = [[  12.065165,  0.000000,  0.000000],
      [  0.000000,  12.065165,  0.000000],
      [  0.000000,  0.000000,  12.065165]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 1
sc_ind_i = 8
sc_ind_j = 9
fc = [[ -5.806100,  0.000000,  0.000000],
      [  0.000000,  0.402342,  0.000000],
      [  0.000000,  0.000000,  0.402342]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 1
sc_ind_i = 8
sc_ind_j = 10
fc = [[  0.402342,  0.000000,  0.000000],
      [  0.000000, -5.806100,  0.000000],
      [  0.000000,  0.000000,  0.402342]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 1
sc_ind_i = 8
sc_ind_j = 11
fc = [[ -0.676651,  0.000000,  0.000000],
      [  0.000000, -0.676651,  0.000000],
      [  0.000000,  0.000000,  0.013278]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 1
sc_ind_i = 8
sc_ind_j = 12
fc = [[  0.402342,  0.000000,  0.000000],
      [  0.000000,  0.402342,  0.000000],
      [  0.000000,  0.000000, -5.806100]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 1
sc_ind_i = 8
sc_ind_j = 13
fc = [[ -0.676651,  0.000000,  0.000000],
      [  0.000000,  0.013278,  0.000000],
      [  0.000000,  0.000000, -0.676651]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 1
sc_ind_i = 8
sc_ind_j = 14
fc = [[  0.013278,  0.000000,  0.000000],
      [  0.000000, -0.676651,  0.000000],
      [  0.000000,  0.000000, -0.676651]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

basis_ind = 1
sc_ind_i = 8
sc_ind_j = 15
fc = [[ -0.407391,  0.000000,  0.000000],
      [  0.000000, -0.407391,  0.000000],
      [  0.000000,  0.000000, -0.407391]]
fc_mats.set_force_constant(basis_ind,sc_ind_i,sc_ind_j,fc,lattice)

print(fc_mats.fc_matrices[0,:,:])

# ---------------------------------------------------------------------------------------------------

"""
Fourier transform the fc matrices to get the dynamical matrices. impose acoustic sum rules etc 
either here or to the fc mats.
"""

dyn_mats = mod_dynmat.dynamical_matrices()

# ---------------------------------------------------------------------------------------------------









