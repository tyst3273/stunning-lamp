import numpy as np
from mod_utils import print_message, raise_error

class lattice:

    """
    hold primitive and supercell vectors, atoms positions, reciprocal lattices, neighbor vectors in
    cartesian and spherical coords, etc.
    """

    # -----------------------------------------------------------------------------------------------

    def __init__(self):

        """
        nothing for now
        """

        pass

    # -----------------------------------------------------------------------------------------------

    # -----------------------------------------------------------------------------------------------
    # public methods go below here
    # -----------------------------------------------------------------------------------------------

    def set_primitive_vectors(self,prim_vecs,prim_scale):

        """
        set primitive lattice vectors and find reciprocal lattice
        """

        # check the primitive vectors
        try:
            self.prim_vecs = np.array(prim_vecs,dtype=float)
        except:
            message = 'prim_vecs should be a 3x3 list of floats'
            raise_error(message)
        if self.prim_vecs.shape[0] !=3 and self.prim_vecs.shape[1] != 3:
            message = 'prim_vecs should be a 3x3 list of floats'
            raise_error(message)

        # check the scale
        try:
            prim_scale = float(prim_scale)
        except:
            message = 'prim_scale should be a positive float'
            raise_error(message)
        if not prim_scale > 0:
            message = 'prim_scale should be a positive float'
            raise_error(message)

        pass

    # -----------------------------------------------------------------------------------------------

    def set_basis_positions(self,basis_pos,basis_labels=None):

        """
        set basis atom positions and convert them to cartesian coordinates
        """

        pass

    # -----------------------------------------------------------------------------------------------

    def set_supercell_vectors(self,sc_vecs,sc_scale):

        """
        set the supercell vectors, find reciprocal lattice, and find transformation matrix to 
        primitive cell
        """

        # check the primitive vectors
        try:
            self.sc_vecs = np.array(sc_vecs,dtype=float)
        except:
            message = 'sc_vecs should be a 3x3 list of floats'
            raise_error(message)
        if self.sc_vecs.shape[0] !=3 and self.sc_vecs.shape[1] != 3:
            message = 'sc_vecs should be a 3x3 list of floats'
            raise_error(message)

        # check the scale
        try:
            sc_scale = float(sc_scale)
        except:
            message = 'sc_scale should be a positive float'
            raise_error(message)
        if not sc_scale > 0:
            message = 'sc_scale should be a positive float'
            raise_error(message)


        pass

    # -----------------------------------------------------------------------------------------------

    def set_supercell_positions(self,sc_pos,sc_labels=None):

        """
        set the supercell atom positions and convert them to cartesian coords.
        """

        pass

    # -----------------------------------------------------------------------------------------------

    def find_neighbors(self):

        """
        loop over atoms in the supercell, find nearest neighbor lists and vectors, convert vectors
        to cartesian and spherical coords. need spherical coords to set up rotation matrices to be
        used later. 
        """

        pass

    # -----------------------------------------------------------------------------------------------

    def write_neighbor_lists(self):

        """
        write the neighbor lists to files PRIM_NEIGHBORS.OUT and SC_NEIGHBORS.OUT PRIM_NN contains 
        the neighbor lists between the basis atoms and sc atoms. SC_NEIGHBORS contains neighbor 
        lists for all atoms in the supercell. 
        """

        pass

    # -----------------------------------------------------------------------------------------------

    # -----------------------------------------------------------------------------------------------
    # private methods go below here
    # -----------------------------------------------------------------------------------------------



