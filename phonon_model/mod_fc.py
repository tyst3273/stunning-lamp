import numpy as np
from mod_utils import raise_error

class force_constant_matrices:
    
    # -------------------------------------------------------------------------------------------------

    def __init__(self,lattice):

        # fc matrices between each atom and the atoms in the SC
        self.fc_matrices = np.zeros((lattice.num_basis,lattice.num_basis,lattice.num_sc,3,3))

        # copy some stuff
        self.num_basis = lattice.num_basis
        self.num_sc = lattice.num_sc

    # -------------------------------------------------------------------------------------------------

    def set_force_constant(self,basis_ind_i,basis_ind_j,sc_ind_i,sc_ind_j,fc,lattice):

        """
        set force constants. there should be 1 between each primitive cell atom and supercell atom.

        fc_ij = [[ phi(i,j)_x,x , phi(i,j)_x,y , phi(i,j)_x,z ],
                 [ phi(i,j)_y,x , phi(i,j)_y,y , phi(i,j)_y,z ],
                 [ phi(i,j)_z,x , phi(i,j)_z,y , phi(i,j)_z,z ]]

        fc units = eV/(Angst*Angstr)
        """

        fc = np.array(fc)
        if fc.shape[0] != 3 or fc.shape[1] != 3:
            message = 'force constant matrix should be 3x3'
            raise_error(message)

        self.fc_matrices[basis_ind_i,basis_ind_j,sc_ind_j,:,:] = np.array(fc)

    # -------------------------------------------------------------------------------------------------

    def parse_force_constants_file(self,fc_file='force_constants'):

        with open(fc_file,'r') as f_fc:
            
            for ii in range(self.num_basis):
                for jj in range(self.num_sc):
                    inds = f_fc.readline().strip().split()
                    inds = [int(x) for x in inds]
                    basis_ind_i = inds[0]
                    basis_ind_j = inds[1]
                    sc_ind_i = inds[2]
                    sc_ind_j = inds[3]
                    fc = np.zeros((3,3))
                    for kk in range(3):
                        line = f_fc.readline().strip().split()
                        fc[kk,:] = [float(x) for x in line]
                    self.fc_matrices[basis_ind_i,basis_ind_j,sc_ind_j,:,:] = fc
                
    # -------------------------------------------------------------------------------------------------




