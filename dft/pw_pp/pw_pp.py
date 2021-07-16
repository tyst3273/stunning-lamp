import mod_crystal
import mod_basis 
import mod_kpoints
#import mod_pp
#import mod_hamiltonian
#import mod_scf

# -----------------------------------------------------------------------------------------------------

crystal = mod_crystal.crystal()

# set lattice vectors and compute reciprocal lattice
lat_vecs = [[0.0,1.0,1.0],      # in angstrom
            [1.0,0.0,1.0],
            [1.0,1.0,0.0]]
lat_scale = 2.01                # scale for lattice vectors
crystal.set_lattice(lat_vecs,lat_scale)

# set atom positions and convert to cartesian coords
atom_pos = [[0.00,0.00,0.00],   # in crystal coords
            [0.25,0.25,0.25]]
atom_types = [1,1]              # 1 per atom
atomic_numbers = [14]           # 1 per type
crystal.set_atoms(atom_pos,atom_types,atomic_numbers)

# -----------------------------------------------------------------------------------------------------

pw = mod_basis.pw()

# set up plane waves, max G vectors, etc
pw_cut_off = 250 # eV
pw.set_plane_waves(crystal,pw_cut_off)

# -----------------------------------------------------------------------------------------------------

kpts = mod_kpoints.kpoints()

# set kpoint grid 
kpt_grid = [7,7,7]
grid_reduction = 'full' 
kpts.set_grid(kpt_grid,grid_reduction)
 
# -----------------------------------------------------------------------------------------------------






