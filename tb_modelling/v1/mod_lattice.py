import numpy as np
from mod_utils import print_stdout, raise_error


class lattice:

    # ------------------------------------------------------------------------------------------------

    def __init__(self):

        pass

    # ------------------------------------------------------------------------------------------------

    def set_atom_data(self,atom_masses=[51.9961],atomic_numbers=[24]):

        """
        set data for atoms. atomic masses, atomic numbers, etc.
        should be lists with 1 entry for each type given in set_*_pos() methods. 
        """
    
        self.atom_masses = atom_masses
        self.atomic_numbers = atomic_numbers

    # ------------------------------------------------------------------------------------------------
        
    # ================================================================================================
    # ----------------------------------- primitive cell ---------------------------------------------
    # ================================================================================================

    def set_prim_vecs(self,prim_vecs=[[1,0,0],[0,1,0],[0,0,1]],prim_scale=2.8525):

        """
        set primitive lattice vectors. prim_vecs should be a 3x3 matrix in Angstrom, 
        prim_scale is just a scaling factor that multiplies prim_vecs
        """

        # set the primitive lattice vectors
        self.prim_vecs = np.array(prim_vecs)
        if self.prim_vecs.shape[0] != 3 or self.prim_vecs.shape[1] != 3:
            message = 'prim_vecs should be a 3x3 list of floats'
            raise_error(message)
        self.prim_vecs = self.prim_vecs*prim_scale

        # get the reciprocal lattice vectors
        self.r_prim_vecs, self.prim_vol = self._compute_reciprocal_lattice(self.prim_vecs)

    # -----------------------------------------------------------------------------------------------

    def set_basis_pos(self,basis_pos=[[0.0,0.0,0.0],[0.5,0.5,0.5]],basis_types=[1,1]):

        # basis pos in crystal coords
        self.basis_pos = np.array(basis_pos)
        if self.basis_pos.shape[1] != 3:
            message = 'basis positions should be 3D'
            raise_error(message)
        self.basis_pos_cart = self._transform_pos(self.prim_vecs,self.basis_pos)

        # check that number of atoms matches number of given types
        self.num_basis = self.basis_pos.shape[0]
        if self.num_basis != len(basis_types):
            message = 'number of positions doenst match number of basis types'
            raise_error(message)
        self.basis_types = basis_types

        # check number of types mathces number of masses and number of atomic numbers
        self.num_types = np.unique(self.basis_types).shape[0]
        if self.num_types != len(self.atom_masses):
            message = 'number of unique basis types doesnt match number of masses'
            raise_error(message)
        if self.num_types != len(self.atomic_numbers):
            message = 'number of unique basis types doesnt match number of atomic numbers'
            raise_error(message)

        # set the masses and atomic numbers for the atoms
        self.basis_masses = np.zeros(self.num_basis)
        for atom in range(self.num_basis):
            self.basis_masses[atom] = self.atom_masses[basis_types[atom]-1]
        self.basis_z = np.zeros(self.num_basis)
        for atom in range(self.num_basis):
            self.basis_z[atom] = self.atomic_numbers[basis_types[atom]-1]

    # -----------------------------------------------------------------------------------------------

    # ================================================================================================
    # -------------------------------------- supercell -----------------------------------------------
    # ================================================================================================

    def set_sc_vecs(self,sc_vecs=[[1,0,0],[0,1,0],[0,0,1]],sc_scale=5.7049):

        """
        set supercell lattice vectors. sc_vecs should be a 3x3 matrix in Angstrom,
        sc_scale is just a scaling factor that multiplies sc_vecs
        """

        # set the supercell lattice vectors
        self.sc_vecs = np.array(sc_vecs)
        if self.sc_vecs.shape[0] != 3 or self.sc_vecs.shape[1] != 3:
            message = 'sc_vecs should be a 3x3 list of floats'
            raise_error(message)
        self.sc_vecs = self.sc_vecs*sc_scale

        # get the supercell reciprocal lattice vectors
        self.r_sc_vecs, self.sc_vol = self._compute_reciprocal_lattice(self.sc_vecs)

    # -----------------------------------------------------------------------------------------------

    def set_sc_pos(self,sc_pos=[[0.0,0.0,0.0],[0.5,0.5,0.5]],sc_types=[1,1]):

        # supercell pos in crystal coords
        self.sc_pos = np.array(sc_pos)
        if self.sc_pos.shape[1] != 3:
            message = 'supercell positions should be 3D'
            raise_error(message)
        self.sc_pos_cart = self._transform_pos(self.sc_vecs,self.sc_pos)

        # check that number of atoms matches number of given type
        self.num_sc = self.sc_pos.shape[0]
        if self.num_sc != len(sc_types):
            message = 'number of sc positions doenst match number of sc types'
            raise_error(message)
        self.sc_types = sc_types

        # check number of types mathces number of types in basis 
        num_types = np.unique(self.sc_types).shape[0]
        if num_types != self.num_types:
            message = 'number of unique sc types doesnt match number of basis types'
            raise_error(message)

        # set the masses and atomic numbers for the atoms
        self.sc_masses = np.zeros(self.num_sc)
        for atom in range(self.num_sc):
            self.sc_masses[atom] = self.atom_masses[sc_types[atom]-1]
        self.sc_z = np.zeros(self.num_sc)
        for atom in range(self.num_sc):
            self.sc_z[atom] = self.atomic_numbers[sc_types[atom]-1]

    # -----------------------------------------------------------------------------------------------

    def find_nearest_neighbors(self,basis_inds=None):

        """
        get relative position vectors (using minimum image) between basis atoms and supercell atoms
        and convert those to cartesian and spherical coords
        """

        # indices of basis atoms in the supercell positions
        if basis_inds == None:
            self.basis_inds = np.arange(self.num_basis)
        else:
            self.basis_inds = np.array(basis_inds)
            if self.basis_inds.min() < 0 or self.basis_inds.max() >= self.num_sc:
                message = 'basis_inds should be integers in range [0,number of atoms in supercell]'
                raise_error(message)

        # relative pos vectors. same vector, different coordinate basis
        self.sc_rel = np.zeros((self.num_basis,self.num_sc,3)) # crystal coords
        self.sc_rel_cart = np.zeros((self.num_basis,self.num_sc,3)) # cartesian coords
        self.sc_rel_sph = np.zeros((self.num_basis,self.num_sc,3)) # spherical coords

        message = 'check sign of relative position vectors'
        print_stdout(message,msg_type='WARNING')

        for ii in range(self.num_basis):

            # index of basis atom
            atom = self.basis_inds[ii]

            pos = self.sc_pos[atom,:] # position of this basis atom in the supercell crystal coords
            rel_pos = np.copy(self.sc_pos) # relative positions
            rel_pos[:,0] = rel_pos[:,0]-pos[0]
            rel_pos[:,1] = rel_pos[:,1]-pos[1]
            rel_pos[:,2] = rel_pos[:,2]-pos[2]
            shift = np.zeros((self.num_sc,3))

            # check whether to shift atoms in the a1 direction
            shift[:,0] = -(rel_pos[:,0] >= 1/2).astype(int)
            shift[:,0] = shift[:,0]+(rel_pos[:,0] < -1/2).astype(int)

            # check whether to shift atoms in the a2 direction
            shift[:,1] = -(rel_pos[:,1] >= 1/2).astype(int)
            shift[:,1] = shift[:,1]+(rel_pos[:,1] < -1/2).astype(int)

            # check whether to shift atoms in the a3 direction
            shift[:,2] = -(rel_pos[:,2] >= 1/2).astype(int)
            shift[:,2] = shift[:,2]+(rel_pos[:,2] < -1/2).astype(int)

            # shift relative positions to minimum image
            rel_pos = rel_pos+shift
            self.sc_rel[ii,:,:] = rel_pos

            # convert relative position in crystal coords to cartesian
            self.sc_rel_cart[ii,:,:] = self._transform_pos(self.sc_vecs,self.sc_rel[ii,:,:])

            # convert relative position from cartesian to spherical coords
            self.sc_rel_sph[ii,:,:] = self._cart_to_spherical(self.sc_rel_cart[ii,:,:])

    # -----------------------------------------------------------------------------------------------

    def _cart_to_spherical(self,cart_pos):

        """
        convert all vectors in the [n_atom x 3] array from cartesian to spherical coords
        """

        cart_pos = np.copy(cart_pos) # need a copy, not pointer. it was misbehaving
        sph_coords = np.zeros((cart_pos.shape[0],3)) # r, theta, phi

        r = np.sqrt(np.sum(cart_pos**2,axis=1)) # radial distance
        sph_coords[:,0] = np.copy(r)

        # need to handle r == 0 correctly. set arg == 1 so that theta = 0
        ind = np.argwhere(r == 0).flatten()
        r[ind] = 1
        cart_pos[ind,2]
        sph_coords[:,1] = np.arccos(cart_pos[:,2]/r)  # theta

        # atan is defined for arg -> inf
        sph_coords[:,2] = np.arctan2(cart_pos[:,1],cart_pos[:,0])   # phi

        return sph_coords

    # -----------------------------------------------------------------------------------------------
    
    # ===============================================================================================
    # ------------------------------------ helper methods -------------------------------------------
    # ===============================================================================================

    def _transform_pos(self,t_mat,old_pos):

        """
        wrapper to just matrix multiply each atoms vector by the fiven transformation matrix.
        """

        nat = old_pos.shape[0]
        new_pos = np.zeros((nat,3))
        for atom in range(nat):
            new_pos[atom,:] = np.matmul(t_mat,old_pos[atom,:])

        return new_pos

    # -----------------------------------------------------------------------------------------------

    def _compute_reciprocal_lattice(self,lat_vecs):

        """
        convert lattice vectors to reciprocal lattice vectors and return r_lat_vecs and cell volume
        """

        r_lat_vecs = np.zeros((3,3))
        cell_vol = lat_vecs[0,:].dot(np.cross(lat_vecs[1,:],lat_vecs[2,:]))
        r_lat_vecs[0,:] = 2*np.pi*np.cross(lat_vecs[1,:],lat_vecs[2,:])/cell_vol
        r_lat_vecs[1,:] = 2*np.pi*np.cross(lat_vecs[2,:],lat_vecs[0,:])/cell_vol
        r_lat_vecs[2,:] = 2*np.pi*np.cross(lat_vecs[0,:],lat_vecs[1,:])/cell_vol

        return r_lat_vecs, cell_vol

    # -----------------------------------------------------------------------------------------------

    def print_lattice(self):

        """
        print all of the stuff setup in this class to screen for diagnostics
        """

        # print primitive cell
        message = 'real space primitive cell (Angstrom):\n'
        for ii in range(3):
            message = message+f'  {self.prim_vecs[ii,0]: 2.4f} {self.prim_vecs[ii,1]: 2.4f}' \
                    f' {self.prim_vecs[ii,2]: 2.4f}\n'
        message = message+f' primitive cell volume (Angstrom**3):\n   {self.prim_vol:2.3f}\n'
        message = message+' primtive cell reciprocal lattice (1/Angstrom):\n'
        for ii in range(3):
            message = message+f'  {self.r_prim_vecs[ii,0]: 2.4f} {self.r_prim_vecs[ii,1]: 2.4f}' \
                    f' {self.r_prim_vecs[ii,2]: 2.4f}\n'
        message = message+' primitive cell atom positions (crystal coords):\n'
        message = message+'   ind, type,   x  y  z\n'
        for ii in range(self.num_basis):
            message = message+f'   {ii}, {self.basis_types[ii]},  {self.basis_pos[ii,0]: 2.3f}' \
                    f' {self.basis_pos[ii,1]: 2.3f} {self.basis_pos[ii,2]: 2.3f}\n'
        message = message+' primitive cell atom positions (cartesian coords, Angstrom):\n'
        message = message+'   ind, type,   x  y  z\n'
        for ii in range(self.num_basis):
            message = message+f'   {ii}, {self.basis_types[ii]},  {self.basis_pos_cart[ii,0]: 2.3f}' \
                    f' {self.basis_pos_cart[ii,1]: 2.3f} {self.basis_pos_cart[ii,2]: 2.3f}\n'
        print_stdout(message,msg_type='PRIMITIVE CELL')

        # print supercell
        message = 'real space supercell (Angstrom):\n'
        for ii in range(3):
            message = message+f'  {self.sc_vecs[ii,0]: 2.4f} {self.sc_vecs[ii,1]: 2.4f}' \
                    f' {self.sc_vecs[ii,2]: 2.4f}\n'
        message = message+f' supercell volume (Angstrom**3):\n   {self.sc_vol:2.3f}\n'
        message = message+' supercell reciprocal lattice (1/Angstrom):\n'
        for ii in range(3):
            message = message+f'  {self.r_sc_vecs[ii,0]: 2.4f} {self.r_sc_vecs[ii,1]: 2.4f}' \
                    f' {self.r_sc_vecs[ii,2]: 2.4f}\n'
        message = message+' supercell atom positions (crystal coords):\n'
        message = message+'   ind, type,   x  y  z\n'
        for ii in range(self.num_sc):
            message = message+f'   {ii}, {self.sc_types[ii]},  {self.sc_pos[ii,0]: 2.3f}' \
                    f' {self.sc_pos[ii,1]: 2.3f} {self.sc_pos[ii,2]: 2.3f}\n'
        message = message+' supercell atom positions (cartesian coords, Angstrom):\n'
        message = message+'   ind, type,   x  y  z\n'
        for ii in range(self.num_sc):
            message = message+f'   {ii}, {self.sc_types[ii]},  {self.sc_pos_cart[ii,0]: 2.3f}' \
                    f' {self.sc_pos_cart[ii,1]: 2.3f} {self.sc_pos_cart[ii,2]: 2.3f}\n'
        print_stdout(message,msg_type='SUPERCELL')

    # -----------------------------------------------------------------------------------------------

    def print_neighbors(self):

        """
        print neighbor vectors to screen for diagnostics
        """

        # print relative pos
        message = f' there are {self.num_basis} atoms in the primitive cell\n'
        for ii in range(self.num_basis):
            atom = self.basis_inds[ii]
            message = message+f'  basis atom index: {atom},  type: {self.sc_types[atom]},'\
                    f' pos (Angstrom): {self.sc_pos_cart[atom,0]: 2.3f}'\
                    f' {self.sc_pos_cart[atom,1]: 2.3f} {self.sc_pos_cart[atom,2]: 2.3f}\n'
            message = message+'  vector from basis atom to atoms in supercell (crystal coords):\n'
            message = message+f'   ind, type,   x  y  z\n'
            for jj in range(self.num_sc):
                message = message+f'   {jj}, {self.sc_types[jj]},  {self.sc_rel[ii,jj,0]: 2.3f}' \
                    f' {self.sc_rel[ii,jj,1]: 2.3f} {self.sc_rel[ii,jj,2]: 2.3f}\n'
            message = message+'  vector from basis atom to atoms in supercell (cartesian coords):\n'
            message = message+f'   ind, type,   x  y  z\n'
            for jj in range(self.num_sc):
                message = message+f'   {jj}, {self.sc_types[jj]},  {self.sc_rel_cart[ii,jj,0]: 2.3f}' \
                    f' {self.sc_rel_cart[ii,jj,1]: 2.3f} {self.sc_rel_cart[ii,jj,2]: 2.3f}\n'
            message = message+'  vector from basis atom to atoms in supercell (spherical coords):\n'
            message = message+f'   ind, type,   r  theta  phi\n'
            for jj in range(self.num_sc):
                message = message+f'   {jj}, {self.sc_types[jj]},  {self.sc_rel_sph[ii,jj,0]: 2.3f}' \
                    f' {self.sc_rel_sph[ii,jj,1]: 2.3f} {self.sc_rel_sph[ii,jj,2]: 2.3f}\n'
        print_stdout(message,msg_type='NEIGHBORS')

    # -----------------------------------------------------------------------------------------------






