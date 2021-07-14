import mod_lattice
import mod_fc
import mod_dynmat 

# ---------------------------------------------------------------------------------------------------

"""

set up the lattice. need primitive cell positions etc to set up crystal lattice. need supercell to set
up coupling to neighbors

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





