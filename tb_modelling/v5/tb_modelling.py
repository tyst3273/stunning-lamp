# system modules
import sys

# custom modules
import mod_lattice
import mod_atoms
from mod_utils import raise_error

# -------------------------------------------------------------------------------------------------
# get which task to do

if len(sys.argv) != 1:
    try:
        task = int(sys.argv[1])
    except:
        message = f'task \'{task}\' is not one of the allowed options'
        raise_error(task)
else:
    task = 1 # do band structure calculation

# -------------------------------------------------------------------------------------------------
# set up primitive cell

lattice = mod_lattice.lattice()

prim_scale = 3.929 # scaling factor for primitive vectors
prim_vecs = [[1.000, 0.000, 0.000], # row vectors in angstroms
             [0.000, 1.000, 0.000],
             [0.000, 0.000, 2.496]]
lattice.set_prim_cell(prim_vecs,prim_scale)

basis_atoms = []
basis_inds = []
basis_type = []

sc_atoms = []
sc_inds = []
sc_types = []

# -------------------------------------------------------------------------------------------------
# set atoms in primitive cell

atoms = mod_atoms.atoms(lattice)

reduced_pos = [0.0, 0.0, 0.0]
atom_type = 'Cu'

atoms.set_atom(reduced_pos,atom_type)






