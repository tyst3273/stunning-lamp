# system modules
import numpy as np

# custom modules
from mod_utils import raise_error

# -------------------------------------------------------------------------------------------------

class lattice:

    """
    stores primitive and supercell lattice vectors, converts atom positions to cartesian coords, 
    etc. also computes neighbor vectors using minimum image convention
    """

    # ---------------------------------------------------------------------------------------------

    def __init__(self):

        """
        nothing yet
        """

        pass

    # ---------------------------------------------------------------------------------------------

    def set_prim_cell(self,prim_vecs,prim_scale):

        """
        set the primitive cell lattice vectors and get the reciprocal lattice
        """

        try:
            self.prim_vecs = np.array(prim_vecs)*prim_scale
        except:
            message = 'primitive vectors should be 3x3 list of floats'
            raise_error(message)
        if self.prim_vecs.shape[0] != 3 and self.prim_vecs.shape[1] != 3:
            message = 'primitive vectors should be 3x3 list of floats'
            raise_error(message)

        self.r_prim_vecs, self.prim_cell_vol = self._compute_reciprocal_lattice(self.prim_vecs)

    # ---------------------------------------------------------------------------------------------

    def convert_to_cart(self,reduced_position):

        return np.matmul()

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





