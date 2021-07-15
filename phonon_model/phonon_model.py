import mod_lattice
import mod_fc
import mod_qpoints
import mod_dynmat 

# ---------------------------------------------------------------------------------------------------

lattice = mod_lattice.lattice()

### set atom data
atom_masses = [51.9961] 
lattice.set_atom_data(atom_masses)

### read primtive cell vectors and positions
prim_cell_file = 'prim_cell'
lattice.set_primitive_cell(prim_cell_file)

### read super cell vectors and positions
super_cell_file = 'super_cell'
lattice.set_super_cell(super_cell_file)

### find nearest neighbor vectors 
lattice.find_neighbors()

### write lattice info to files
#lattice.write_lattice()
#lattice.write_neighbors()


# -------------------------------------------------------------------------------------------------

fc_mats = mod_fc.force_constant_matrices(lattice)

### read the force constant file 
fc_mats.parse_force_constants_file('force_constants')


# ---------------------------------------------------------------------------------------------------

qpoints = mod_qpoints.qpoints()

### set up q points path
q_vert = [[0.000,0.000,0.000],
          [0.500,0.000,0.000],
          [0.500,0.500,0.000],
          [0.000,0.000,0.000],
          [0.500,0.500,0.500]]
qpoints.q_path(lattice,q_vert,num_q_sm=151)

# ---------------------------------------------------------------------------------------------------

dyn_mats = mod_dynmat.dynamical_matrices(lattice,qpoints)

### diagonalize dynamical mats
dyn_mats.calculate_phonons(lattice,fc_mats,qpoints)

# ---------------------------------------------------------------------------------------------------




