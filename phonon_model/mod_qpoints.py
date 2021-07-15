import numpy as np
from mod_utils import raise_error


class qpoints:

    # ------------------------------------------------------------------------------------------------

    def __init__(self):
    
        pass

    # ------------------------------------------------------------------------------------------------

    def q_path(self,lattice,q_vert=[[0,0,0],[0.5,0,0]],num_q_sm=101):

        # vertices of the q path
        q_vert = np.array(q_vert)
        if q_vert.shape[1] != 3:
            message = 'q vertices should be 3D'
            raise_error(message)
        num_vert = q_vert.shape[0]
        num_path = num_vert-1

        # vertices of the q path in cartesian coords
        q_vert_cart = self._crystal_to_cart(lattice.r_prim_vecs,q_vert)

        self.q_dist = np.zeros(num_path)
        for qq in range(num_path):
            q_dist = q_vert_cart[qq+1,:]-q_vert_cart[qq,:]
            q_dist = np.sqrt((q_dist**2).sum())
            self.q_dist[qq] = q_dist

        num_per = self.q_dist/self.q_dist.min()*num_q_sm
        num_per = np.round(num_per).astype(int)
        
        self.num_q = num_per.sum()
        self.q_cart = np.zeros((self.num_q,3))

        shift = 0
        for qq in range(num_path):
            for ii in range(3):
                self.q_cart[shift:shift+num_per[qq],ii] = np.linspace(q_vert_cart[qq,ii],
                        q_vert_cart[qq+1,ii],num_per[qq])
            shift = shift+num_per[qq]

    # -----------------------------------------------------------------------------------------------

    def _crystal_to_cart(self,lat_vecs,crystal_coords):

        """
        wrapper to just matrix multiply each atoms vector by the fiven transformation matrix.
        """

        nat = crystal_coords.shape[0]
        cart_coords = np.zeros((nat,3))
        for atom in range(nat):
            cart_coords[atom,:] = np.matmul(lat_vecs,crystal_coords[atom,:])

        return cart_coords

    # ------------------------------------------------------------------------------------------------



