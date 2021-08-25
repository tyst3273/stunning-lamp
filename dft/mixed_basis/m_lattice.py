# system modules
import numpy as np

# custom modules
from m_utils import raise_error

class lattice:

    def __init__(self):

        pass

    def set_lattice_vectors(self,lat_vecs,lat_scale=1):

        # check that the lattice vectors make sense
        err_msg = 'lat_vecs should be a 3x3 list of floats'
        try:
            self.lat_vecs = np.array(lat_vecs,dtype=float)
        except:
            raise_error(err_msg)
        if len(self.lat_vecs.shape) != 2:
            raise_error(err_msg)
        if self.lat_vecs.shape[0] != 3 or self.lat_vecs.shape[1] != 3:
            raise_error(err_msg)

        # check that the lat_scale makes sense
        err_msg = 'lat_scale should be a positive, non-zero float'
        try:
            lat_scale = float(lat_scale)
        except:
            raise_error(err_msg)
        if lat_scale <= 0: 
            raise_error(err_msg)

        # now apply the scale
        self.lat_vecs = self.lat_vecs*lat_scale

        





