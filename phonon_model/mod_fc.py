import numpy as np
from mod_utils import raise_error

class force_constant_matrices:
    
    # -------------------------------------------------------------------------------------------------

    def __init__(self,lattice):

        # fc matrices between each atom and the atoms in the SC
        self.fc_matrices = np.zeros((lattice.num_basis*lattice.num_sc,3,3))
        self.inds = np.zeros((lattice.num_basis*lattice.num_sc,4),dtype=int)

        # copy some stuff
        self.num_basis = lattice.num_basis
        self.num_sc = lattice.num_sc

    # -------------------------------------------------------------------------------------------------

    def parse_force_constants_file(self,fc_file='force_constants'):

        with open(fc_file,'r') as f_fc:
            
            for pair in range(self.num_sc*self.num_basis):
                    
                # get the inds from the file
                try:
                    inds = f_fc.readline().strip().split() 
                    inds = [int(x) for x in inds]
                except:
                    message = f'indices in file \'{fc_file}\ seem wrong'
                    raise_error(message)
                basis_ind_i = inds[0]
                basis_ind_j = inds[1]
                sc_ind_i = inds[2]
                sc_ind_j = inds[3]

                # get the fc from the file
                fc = np.zeros((3,3))
                try:
                    for kk in range(3):
                        line = f_fc.readline().strip().split()
                        fc[kk,:] = [float(x) for x in line]
                except:
                    message = 'force constants in file \'{fc_file}\ seem wrong'
                    raise_error(message)

                self.fc_matrices[pair,:,:] = fc
                self.inds[pair,:] = [basis_ind_i,basis_ind_j,sc_ind_i,sc_ind_j]

    # -------------------------------------------------------------------------------------------------




