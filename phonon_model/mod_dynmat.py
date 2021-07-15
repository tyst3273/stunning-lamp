import numpy as np
from mod_utils import raise_error

class dynamical_matrices:

    # --------------------------------------------------------------------------------------------------

    def __init__(self,lattice,qpoints):

        # dynamical matrices
        self.dyn_mats = np.zeros((qpoints.num_q,lattice.num_basis*3,lattice.num_basis*3),dtype=complex)

        # eigenvals
        self.freq = np.zeros((qpoints.num_q,lattice.num_basis*3))
        self.eigenvecs = np.zeros((qpoints.num_q,lattice.num_basis*3,lattice.num_basis*3),dtype=complex)

    # --------------------------------------------------------------------------------------------------

    def calculate_phonons(self,lattice,fc_mats,qpoints):

        for qq in range(qpoints.num_q):

            q_vec = qpoints.q_cart[qq,:]
            self.dyn_mats[qq,:,:] = self._compute_dynamical_matrix(q_vec,fc_mats,lattice)

            self.freq[qq,:] = np.sqrt(np.abs(np.linalg.eigvalsh(self.dyn_mats[qq,:,:])))

        np.savetxt('freqs',self.freq,fmt='%2.4f')

    # --------------------------------------------------------------------------------------------------

    def _compute_dynamical_matrix(self,q_vec,fc_mats,lattice):

        dyn_mat = np.zeros((lattice.num_basis*3,lattice.num_basis*3),dtype=complex)

        for pair in range(lattice.num_sc*lattice.num_basis):

            ii = fc_mats.inds[pair,0]
            jj = fc_mats.inds[pair,1]
            type_i = lattice.basis_types[ii]
            type_j = lattice.basis_types[jj]
            sc_ii = fc_mats.inds[pair,2]
            sc_jj = fc_mats.inds[pair,3]

            fc = fc_mats.fc_matrices[pair,:,:]
            q_dot_R = lattice.sc_rel_cart[sc_ii,sc_jj]*q_vec
            mass_fac = 1/np.sqrt(lattice.atom_masses[type_i]*lattice.atom_masses[type_j])
            fc = fc*np.exp(1j*q_dot_R.sum())/mass_fac
            dyn_mat[ii*3:(ii+1)*3,jj*3:(jj+1)*3] = dyn_mat[ii*3:(ii+1)*3,jj*3:(jj+1)*3]+fc

        return dyn_mat

    # --------------------------------------------------------------------------------------------------


