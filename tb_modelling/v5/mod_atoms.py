# system modules
import numpy as np

# custom modules
from mod_utils import raise_error

# -------------------------------------------------------------------------------------------------

class _atom:
    
    """
    atom object. it holds each atom's positions, orbitals, etc
    """

    def __init__(self):

        pass

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

class atoms:

    """
    set atom positions, orbitals, etc
    """

    # ---------------------------------------------------------------------------------------------

    def __init__(self,lattice):

        """
        hold list of 'atoms' and number of atoms
        """

        self.lattice = lattice 

        self.basis_atoms = []
        self.num_basis_atoms = 0
        self.basis_types = []
        self.num_basis_types = 0

        self.sc_atoms = []
        self.num_sc_atoms = 0
        self.sc_types = []
        self.num_sc_types = 0

        self.num_atoms = 0

    # ---------------------------------------------------------------------------------------------

    def set_atom(self,reduced_pos,atom_type='none'):

        """
        creates an instance of '_atom' class and sets some attributes
        """

        try:
            self.reduced_pos = np.array(reduced_pos).flatten()
        except:
            message = 'reduced position should be a 3 element list of floats'
            raise_error(message)
        if self.reduced_pos.size != 3:
            message = 'reduced position should be a 3 element list of floats'
            raise_error(message)

        atom = _atom()
        atom.reduced_pos = reduced_pos
        atom.atom_type = atom_type
        
        if self.reduced_pos[0] >= 1 or self.reduced_pos[1] >= 1 or self.reduced_pos[2] >= 1:
            atom.is_basis = True
            atom.atom_index = self.num_atoms
            self.basis_atoms.append(atom)
            self.num_atoms = self.num_atoms+1
            self.num_basis_atoms = self.num_basis_atoms+1
        else:
            atom.is_basis = False
            atom.atom_index = self.num_atoms
            self.sc_atoms.append(atom)
            self.num_atoms = self.num_atoms+1
            self.num_sc_atoms = self.num_sc_atoms+1

    # ---------------------------------------------------------------------------------------------


    # ---------------------------------------------------------------------------------------------


    # ---------------------------------------------------------------------------------------------
