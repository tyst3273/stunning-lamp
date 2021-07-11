import numpy as np
import mod_io

# ===================================================================================================
# ---------------------------------------------------------------------------------------------------
# ===================================================================================================

class potential:

    # ---------------------------------------------------------------------------------------------

    def __init__(self,lattice,num_pair_types):

        # copy from sc data for ease of use later
        self.basis_types = np.array(lattice.sc_basis_types)

        # keep count of number of pairs that are defined
        self.num_pair_types = num_pair_types

        # data structure to hold all the pair potential parameters, eg force-constants, hopping elements
        self.pair_atoms = np.zeros((num_pair_types,2),dtype=int)
        self.pair_matrix = np.zeros((self.num_pair_types,lattice.sc_num_atoms,
                                                    lattice.sc_num_atoms),dtype=int)

    # -----------------------------------------------------------------------------------------------

    def set_pair(self,pair_type,atom_types):

        """
        set up a matrix where elements are non-zero only between the two types in the pair.
        note, it is by TYPE, not by atom. the data is saved in a dictionary for easy lookup later
        """

        # store the atom types
        self.pair_atoms[pair_type,0] = atom_types[0]
        self.pair_atoms[pair_type,1] = atom_types[1]

        types_1 = (self.basis_types == atom_types[0])
        types_2 = (self.basis_types == atom_types[1])

        if atom_types[0] == atom_types[1]:
            self.pair_matrix[pair_type,:,:] = np.outer(types_1,types_2)
        else:
            self.pair_matrix[pair_type,:,:] = np.outer(types_1,types_2)+np.outer(types_2,types_1)

    # ----------------------------------------------------------------------------------------------

# ===================================================================================================
# ---------------------------------------------------------------------------------------------------
# ===================================================================================================

class from_file:

    # ----------------------------------------------------------------------------------------------

    def __init__(self,lattice,fc_file):

        fc_mats = np.loadtxt(fc_file)
        fc_mats = fc_mats.reshape(lattice.sc_num_atoms,4,3)
        fc_mats = fc_mats[:,1:,:]
        self.fc_mats = fc_mats

    # ----------------------------------------------------------------------------------------------

    def q_path(self,qmin=[0,0,0],qmax=[0.5,0,0],num_q=251):

        self.num_q = num_q
        self.q_points = np.zeros((num_q,3))
        self.q_points[:,0] = np.linspace(qmin[0],qmax[0],num_q)
        self.q_points[:,1] = np.linspace(qmin[1],qmax[1],num_q)
        self.q_points[:,2] = np.linspace(qmin[2],qmax[2],num_q)
        
    # ----------------------------------------------------------------------------------------------

    def calculate(self,lattice):

        self.modes = np.zeros((self.num_q,3))

        for q in range(self.num_q):

            rel_pos = lattice.sc_rel_cart[0,:,:]

            q_arr = self.q_points[q,:]
            q_arr = np.matmul(lattice.p_recip_lat_vecs,q_arr) # 1/Angstrom coords
            print(q_arr)
            q_arr = np.tile(q_arr,reps=[lattice.sc_num_atoms,1])

            exp_i_qR = q_arr*rel_pos
            exp_i_qR = np.exp(1j*exp_i_qR.sum(axis=1))
            exp_i_qR = np.tile(exp_i_qR.reshape(lattice.sc_num_atoms,1,1),reps=[1,3,3])

            self.dynmat = np.sum(exp_i_qR*self.fc_mats,axis=0)/lattice.p_basis_masses[0] 

            self.modes[q,:] = np.linalg.eigvalsh(self.dynmat)

    # ----------------------------------------------------------------------------------------------

# ===================================================================================================
# ---------------------------------------------------------------------------------------------------
# ===================================================================================================









