import numpy as np
from mod_utils import raise_error

class force_constant_matrices:
    
    # -------------------------------------------------------------------------------------------------

    def __init__(self,lattice):

        # position vector (in cart coords.) between each atom and the atoms in the SC
        self.rel_pos = np.zeros((lattice.num_basis,lattice.num_sc,3))

        # fc matrices between each atom and the atoms in the SC
        self.fc_matrices = np.zeros((lattice.num_basis,lattice.num_sc*3,3))

    # -------------------------------------------------------------------------------------------------

    def set_force_constant(self,basis_ind,sc_ind_i,sc_ind_j,fc,lattice):

        """
        set force constants. there should be 1 between each primitive cell atom and supercell atom.

        fc_ij = [[ phi(i,j)_x,x , phi(i,j)_x,y , phi(i,j)_x,z ],
                 [ phi(i,j)_y,x , phi(i,j)_y,y , phi(i,j)_y,z ],
                 [ phi(i,j)_z,x , phi(i,j)_z,y , phi(i,j)_z,z ]]

        fc units = eV/(Angst*Angstr)
        """

        # for phi(i,j)_a,b; R_ij = R(j)-R(i). may need to check sign of this
        self.rel_pos[basis_ind,sc_ind_j,:] = lattice.sc_pos_cart[sc_ind_j,:]-\
                lattice.sc_pos_cart[sc_ind_i,:]

        fc = np.array(fc)
        if fc.shape[0] != 3 or fc.shape[1] != 3:
            print(sc_ind_i,sc_ind_j)
            message = 'force constant matrix should be 3x3'
            raise_error(message)

        self.fc_matrices[basis_ind,sc_ind_j*3:(sc_ind_j+1)*3,:] = np.array(fc)

    # -------------------------------------------------------------------------------------------------





