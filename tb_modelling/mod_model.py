# system modules
import numpy as np 

# custom modules
from mod_utils import print_message, raise_error


# -------------------------------------------------------------------------------------------------
# lattice class

class lattice:

    """
    holds lattice vectors and methods to do lattice vectory stuff
    """

    #----------------------------------------------------------------------------------------------

    def __init__(self):

        """
        nothing yet
        """

        pass

    # ---------------------------------------------------------------------------------------------

    def set_lat_vectors(self,lat_vecs,lat_scale):

        """
        set the real space lattice vectors and call method to calculate reciprocal lattice vectors
        """

        try:
            self.lat_vecs = np.array(lat_vecs)*lat_scale
        except:
            message = 'lattice vectors should be a 3x3 list of floats'
            raise_error(message)
        if self.lat_vecs.shape[0] != 3 or self.lat_vecs.shape[1] != 3:
            message = 'lattice vectors should be a 3x3 list of floats'
            raise_error(message)

        self.r_lat_vecs, self.lat_cell_vol = self._compute_reciprocal_lattice(self.lat_vecs)

    # --------------------------------------------------------------------------------------------- 
    # private methods

    def _compute_reciprocal_lattice(self,lat_vecs):

        """
        convert real space to reciprocal lattice vectors
        """

        r_lat_vecs = np.zeros((3,3))
        cell_vol = lat_vecs[0,:].dot(np.cross(lat_vecs[1,:],lat_vecs[2,:]))
        r_lat_vecs[0,:] = 2*np.pi*np.cross(lat_vecs[1,:],lat_vecs[2,:])/cell_vol
        r_lat_vecs[1,:] = 2*np.pi*np.cross(lat_vecs[2,:],lat_vecs[0,:])/cell_vol
        r_lat_vecs[2,:] = 2*np.pi*np.cross(lat_vecs[0,:],lat_vecs[1,:])/cell_vol

        return r_lat_vecs, cell_vol

    # --------------------------------------------------------------------------------------------- 

    def _convert_to_cartesian(self,reduced_position):

        """
        convert from rduced (crystal) to cartesian coordinates
        """

        return np.matmul(self.lat_vecs.T,reduced_pos)

    # ---------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------------
# atoms class

class atoms:
    
    """
    holds atoms, positions, vectors between atoms, etc
    """

    # ---------------------------------------------------------------------------------------------

    def __init__(self,reduced_pos,atom_type,atom_index):

        

        self.reduced_pos = reduced_pos
        self.atom_type = atom_type

    # --------------------------------------------------------------------------------------------- 
   

# -------------------------------------------------------------------------------------------------
# the main 'model' class that interacts with lattice, atoms, etc. this builds the hamiltonian
 
class model:

    """
    tight binding model class. holds hamiltonian, etc. should have methods to set up model, get 
    hamiltonian matrix at a given k-point, etc.
    """

    # ---------------------------------------------------------------------------------------------

    def __init__(self):

        """
        initialize the lattice and atoms classes and container for the atoms
        """

        self.lattice = lattice()

        self.basis_atoms = []
        self.num_basis_atoms = 0

    # ---------------------------------------------------------------------------------------------

    def set_lattice(self,lat_vecs,lat_scale):

        """
        interface with the lattice class to set up the lattice
        """

        self.lattice.set_lat_vectors(lat_vecs,lat_scale)

    # ---------------------------------------------------------------------------------------------

    def set_basis_atom(self,reduced_pos,atom_type):

        """
        interface with atoms class to create an atom, set properties, and add it to container
        """

        atom_index = self.num_basis_atoms
        self.num_basis_atoms = self.num_basis_atoms+1
        basis_atom = self.atoms(reduced_pos,atom_type,atom_index)
        self.basis_atoms.append(basis_atom)

    # ---------------------------------------------------------------------------------------------




