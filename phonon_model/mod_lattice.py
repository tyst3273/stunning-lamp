import numpy as np
from mod_utils import print_stdout, raise_error


class lattice:

    # ------------------------------------------------------------------------------------------------

    def __init__(self):

        pass

    # ------------------------------------------------------------------------------------------------

    def set_atom_data(self,atom_masses=[51.9961]):
    
        self.atom_masses = atom_masses

    # ------------------------------------------------------------------------------------------------
        
    # ================================================================================================
    # ----------------------------------- primitive cell ---------------------------------------------
    # ================================================================================================

    def set_prim_vecs(self,prim_vecs=[[1,0,0],[0,1,0],[0,0,1]],prim_scale=2.8525):

        self.prim_vecs = np.array(prim_vecs)
        if self.prim_vecs.shape[0] != 3 or self.prim_vecs.shape[1] != 3:
            message = 'prim_vecs should be a 3x3 list of floats'
            raise_error(message)

        self.prim_vecs = self.prim_vecs*prim_scale
        self.r_prim_vecs, self.prim_vol = self._compute_reciprocal_lattice(self.prim_vecs)

    # -----------------------------------------------------------------------------------------------

    def set_basis_pos(self,basis_pos=[[0.0,0.0,0.0],[0.5,0.5,0.5]],basis_types=[1,1]):

        self.basis_pos = np.array(basis_pos)
        if self.basis_pos.shape[1] != 3:
            message = 'basis positions should be 3D'
            raise_error(message)
        self.basis_pos_cart = self._transform_pos(self.prim_vecs,self.basis_pos)

        self.num_basis = self.basis_pos.shape[0]
        if self.num_basis != len(basis_types):
            message = 'number of positions doenst match number of basis types'
            raise_error(message)
        self.basis_types = basis_types

        self.num_types = np.unique(self.basis_types).shape[0]
        if self.num_types != len(self.atom_masses):
            message = 'number of unique basis types doesnt match number of masses'
            raise_error(message)

        self.basis_masses = np.zeros(self.num_basis)
        for atom in range(self.num_basis):
            self.basis_masses[atom] = self.atom_masses[basis_types[atom]-1]

    # -----------------------------------------------------------------------------------------------

    # ================================================================================================
    # -------------------------------------- supercell -----------------------------------------------
    # ================================================================================================

    def set_sc_vecs(self,sc_vecs=[[1,0,0],[0,1,0],[0,0,1]],sc_scale=5.7049):

        self.sc_vecs = np.array(sc_vecs)
        if self.sc_vecs.shape[0] != 3 or self.sc_vecs.shape[1] != 3:
            message = 'sc_vecs should be a 3x3 list of floats'
            raise_error(message)

        self.sc_vecs = self.sc_vecs*sc_scale
        self.r_sc_vecs, self.sc_vol = self._compute_reciprocal_lattice(self.sc_vecs)

    # -----------------------------------------------------------------------------------------------

    def set_sc_pos(self,sc_pos=[[0.0,0.0,0.0],[0.5,0.5,0.5]],sc_types=[1,1]):

        self.sc_pos = np.array(sc_pos)
        if self.sc_pos.shape[1] != 3:
            message = 'supercell positions should be 3D'
            raise_error(message)
        self.sc_pos_cart = self._transform_pos(self.sc_vecs,self.sc_pos)

        self.num_sc = self.sc_pos.shape[0]
        if self.num_sc != len(sc_types):
            message = 'number of sc positions doenst match number of sc types'
            raise_error(message)
        self.sc_types = sc_types

        num_types = np.unique(self.sc_types).shape[0]
        if num_types != len(self.atom_masses):
            message = 'number of unique sc types doesnt match number of masses'
            raise_error(message)

        self.sc_masses = np.zeros(self.num_sc)
        for atom in range(self.num_sc):
            self.sc_masses[atom] = self.atom_masses[sc_types[atom]-1]

    # -----------------------------------------------------------------------------------------------

    # ===============================================================================================
    # -------------------------------- methods for lattice ------------------------------------------
    # ===============================================================================================

    def _transform_pos(self,transf_mat,old_pos):

        nat = old_pos.shape[0]
        new_pos = np.zeros((nat,3))
        for atom in range(nat):
            new_pos[atom,:] = np.matmul(transf_mat,old_pos[atom,:])

        return new_pos

    # -----------------------------------------------------------------------------------------------

    def _compute_reciprocal_lattice(self,lat_vecs):

        r_lat_vecs = np.zeros((3,3))
        cell_vol = lat_vecs[0,:].dot(np.cross(lat_vecs[1,:],lat_vecs[2,:]))
        r_lat_vecs[0,:] = 2*np.pi*np.cross(lat_vecs[1,:],lat_vecs[2,:])/cell_vol
        r_lat_vecs[1,:] = 2*np.pi*np.cross(lat_vecs[2,:],lat_vecs[0,:])/cell_vol
        r_lat_vecs[2,:] = 2*np.pi*np.cross(lat_vecs[0,:],lat_vecs[1,:])/cell_vol

        return r_lat_vecs, cell_vol

    # -----------------------------------------------------------------------------------------------









