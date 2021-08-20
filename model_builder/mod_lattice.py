# system modules
import numpy as np

# custom modules
from mod_utils import raise_error

# --------------------------------------------------------------------------------------------------
# the lattice class

class lattice:
    
    """
    lattice and supercell vectors, reciprocal lattice, basis atom positions, etc.
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self):

        """
        nothing yet
        """

        pass

    # ----------------------------------------------------------------------------------------------

    def set_prim_vecs(self,prim_vecs,prim_scale=1):
        
        """
        set up primitive vectors and reciprocal lattice
        """

        try:
            self.prim_vecs = np.array(prim_vecs)*prim_scale
        except:
            message = 'prim_vecs should be a 3x3 list of floats, prim_scale should be a float'
            raise_error(message)
        if len(self.prim_vecs.shape) != 2:
            message = 'prim_vecs should be a 3x3 list of floats'
            raise_error(message)
        if self.prim_vecs.shape[0] != 3 and self.prim_vecs.shape[1] != 3:
            message = 'prim_vecs should be a 3x3 list of floats'
            raise_error(message)
        
        self.r_lat_vecs, self.prim_cell_vol = self._compute_reciprocal_lattice(self.prim_vecs)

    # ----------------------------------------------------------------------------------------------

    def set_basis_positions(self,basis_pos,basis_labels):

        """
        set atom positions, convert to cartesian coords, and find unique basis labels.
        """

        try:
            self.basis_pos = np.array(basis_pos)
        except:
            message = 'basis_pos should be a n_atom x 3 list of floats'
            raise_error(message)
        if len(self.basis_pos.shape) != 2:
            message = 'basis_pos should be a n_atom x 3 list of floats'
            raise_error(message)
        if self.basis_pos.shape[1] != 3:
            message = 'basis_pos should be a n_atom x 3 list of floats'
            raise_error(message)
        self.num_atoms = self.basis_pos.shape[0]

        if len(basis_labels) != self.num_atoms:
            message = 'number of basis labels does not match number of atoms'
            raise_error(message)

        self.basis_cart = self._convert_to_cartesian(self.basis_pos,self.prim_vecs)

        self.basis_labels = basis_labels
        self.unique_labels = []
        self.num_types = 0
        for label in basis_labels:
            if label not in self.unique_labels:
                self.num_types = self.num_types+1
                self.unique_labels.append(label)

    # ----------------------------------------------------------------------------------------------

    def set_supercell(self,sc_vecs,sc_scale=1):

        """
        set up supercell. given supercell vectors, check that the requested supercell can be made
        then set it up, copy atoms, etc.
        """

        try:
            self.sc_vecs = np.array(sc_vecs)*sc_scale
        except:
            message = 'sc_vecs should be a 3x3 list of floats, sc_scale should be a float'
            raise_error(message)
        if len(self.sc_vecs.shape) != 2:
            message = 'sc_vecs should be a 3x3 list of floats'
            raise_error(message)
        if self.sc_vecs.shape[0] != 3 and self.sc_vecs.shape[1] != 3:
            message = 'sc_vecs should be a 3x3 list of floats'
            raise_error(message)

        message = 'lattice.set_supercell(*) is not complete'
        raise_error(message)

    # ----------------------------------------------------------------------------------------------
    # private methods

    def _compute_reciprocal_lattice(self,lat_vecs):

        """
        convert lattice vectors to reciprocal lattice vectors and return r_lat_vecs and cell volume
        """

        r_lat_vecs = np.zeros((3,3))
        cell_vol = lat_vecs[0,:].dot(np.cross(lat_vecs[1,:],lat_vecs[2,:]))
        r_lat_vecs[0,:] = 2*np.pi*np.cross(lat_vecs[1,:],lat_vecs[2,:])/cell_vol
        r_lat_vecs[1,:] = 2*np.pi*np.cross(lat_vecs[2,:],lat_vecs[0,:])/cell_vol
        r_lat_vecs[2,:] = 2*np.pi*np.cross(lat_vecs[0,:],lat_vecs[1,:])/cell_vol

        return r_lat_vecs, cell_vol

    # -----------------------------------------------------------------------------------------------

    def _convert_to_cartesian(self,reduced_pos,lat_vecs):

        num_atoms = reduced_pos.shape[0]
        cart_pos = np.zeros((num_atoms,3))
        for ii in range(num_atoms):
            cart_pos[ii,:] = np.matmul(lat_vecs.T,reduced_pos[ii,:])
        
        return cart_pos

    # -----------------------------------------------------------------------------------------------















