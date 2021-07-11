import numpy as np
import mod_io

# ===================================================================================================
# ---------------------------------------------------------------------------------------------------
# ===================================================================================================

class lattice:

    # -----------------------------------------------------------------------------------------------

    def __init__(self,p_lat_vecs=[[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]],
                p_basis_pos=[[0.0,0.0,0.0],[0.5,0.5,0.5]],p_basis_types=[1,2],
                p_basis_masses=[26.9815,26.9815],sc_matrix=[[1,0,0],[0,1,0],[0,0,1]],
                sc_lat_vecs=None):
    
        self.sc_matrix = sc_matrix
        self.sc_lat_vecs = sc_lat_vecs

        self.p_lat_vecs = np.array(p_lat_vecs)
        self.p_basis_pos = np.array(p_basis_pos)
        self.p_basis_types = np.array(p_basis_types)
        self.p_basis_masses = np.array(p_basis_masses)
        self.num_types = self.p_basis_types.shape[0]

        self._prep_primitive_cell()
        self._prep_supercell()

    # ----------------------------------------------------------------------------------------------

    def _prep_primitive_cell(self):

        self.p_num_atoms = self.p_basis_pos.shape[0]

        # get reciprocal lattice
        self.p_recip_lat_vecs, self.p_cell_vol = self._compute_reciprocal_lattice(self.p_lat_vecs)

        # convert reduced coords to cartesian coords
        self.p_basis_cart = self._crystal_to_cart(self.p_lat_vecs,self.p_basis_pos)

        # get inverse of a
        self.p_lat_inv = np.linalg.inv(self.p_lat_vecs)


    # ------------------------------------------------------------------------------------------------

    def _prep_supercell(self):

        """
        supercell lattice vectors A are determined from a*T where a are the primitive lattice
        vectors and T is the transformation.
        """

        # this is the matrix T above
        self.sc_matrix = np.array(self.sc_matrix)

        # if sc_lat_vecs are given, use those. otherwise calculate from sc_matrix
        if self.sc_lat_vecs == None:
            self.sc_lat_vecs = np.matmul(self.p_lat_vecs,self.sc_matrix) # A = a*T
        else:
            self.sc_lat_vecs = np.array(self.sc_lat_vecs)                     
            self.sc_matrix = np.matmul(self.p_lat_inv,self.sc_lat_vecs)  # a^(-1)*A = T

        # invert A and T
        self.sc_matrix_inv = np.linalg.inv(self.sc_matrix)
        self.sc_lat_inv = np.linalg.inv(self.sc_lat_vecs)

        # get reciprocal lattice
        self.sc_recip_lat_vecs, self.sc_cell_vol = self._compute_reciprocal_lattice(self.sc_lat_vecs)

        # multiplicity of supercell
        self.num_p_per_sc = int(round(np.linalg.det(self.sc_matrix))) # det(T)
        self.sc_num_atoms = self.num_p_per_sc*self.p_num_atoms

        # check is sc_matrix is diagonal or not and add supercell atoms
        if np.count_nonzero(self.sc_matrix-np.diag(np.diagonal(self.sc_matrix))) != 0:
            self.sc_diag = False
            self._fill_supercell()
        else:
            self.sc_diag = True
            self._fill_diagonal_supercell()

        # check if supercell lattice vectors are orthorhombic
        if np.count_nonzero(self.sc_lat_vecs-np.diag(np.diagonal(self.sc_lat_vecs))) != 0:
            self.sc_ortho = False
        else:
            self.sc_ortho = True
        
        # convert reduced coords to cartesian coords
        self.sc_basis_cart = self._crystal_to_cart(self.sc_lat_vecs,self.sc_basis_pos) 

        # find relative positions between primitive cell atoms and supercell atoms
        self.sc_rel_pos = np.zeros((self.p_num_atoms,self.sc_num_atoms,3))
        self.sc_rel_cart = np.zeros((self.p_num_atoms,self.sc_num_atoms,3)) # x, y, z
        self.sc_rel_sph = np.zeros((self.p_num_atoms,self.sc_num_atoms,3)) # r, theta, phi
        self._minimum_image()

        print(self.sc_basis_cart)
        print(self.sc_rel_pos)
        print(self.sc_rel_cart)

    # -----------------------------------------------------------------------------------------------

    def _crystal_to_cart(self,lat_vecs,old_pos):

        nat = old_pos.shape[0]
        new_pos = np.zeros((nat,3))
        for atom in range(nat):
            new_pos[atom,:] = np.matmul(lat_vecs,old_pos[atom,:])

        return new_pos

    # -----------------------------------------------------------------------------------------------

    def _cart_to_spherical(self,cart_pos):

        cart_pos = np.copy(cart_pos) # need a copy, not pointer
        sph_coords = np.zeros((cart_pos.shape[0],3)) # r, theta, phi

        r = np.sqrt(np.sum(cart_pos**2,axis=1))       
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

    def _compute_reciprocal_lattice(self,lat_vecs):

        r_lat_vecs = np.zeros((3,3))
        cell_vol = lat_vecs[0,:].dot(np.cross(lat_vecs[1,:],lat_vecs[2,:]))
        r_lat_vecs[0,:] = 2*np.pi*np.cross(lat_vecs[1,:],lat_vecs[2,:])/cell_vol
        r_lat_vecs[1,:] = 2*np.pi*np.cross(lat_vecs[2,:],lat_vecs[0,:])/cell_vol
        r_lat_vecs[2,:] = 2*np.pi*np.cross(lat_vecs[0,:],lat_vecs[1,:])/cell_vol

        return r_lat_vecs, cell_vol

    # -----------------------------------------------------------------------------------------------

    def _fill_supercell(self):

        """
        dont know how to do it (yet) in a way thats not wacky-hacky 
        """
        message = '\n *** ERROR ***\n  dont know how to fill non-diagonal supercells yet\n '
        print(message)

        # these are stubs to write the lattice file
        self.sc_basis_pos = np.zeros((self.sc_num_atoms,3))
        self.sc_basis_types = []
        for atom in range(self.sc_num_atoms):
            self.sc_basis_types.append(1)
    
        # write lattice info to file
        mod_io.write_lattice(self)

        exit() # now crash

    # -----------------------------------------------------------------------------------------------

    def _fill_diagonal_supercell(self):

        """
        it is simple for diagonal sc_matrix. just copy basis N1 times along a1, N2 along a2, etc ...
        """

        # set up superlattice basis
        self.sc_basis_pos = np.zeros((self.sc_num_atoms,3))
        self.sc_basis_types = []
    
        count = 0
        shift = np.zeros((self.p_num_atoms,3))
        for ii in range(self.sc_matrix[0,0]):
            for jj in range(self.sc_matrix[1,1]):
                for kk in range(self.sc_matrix[2,2]):
                    shift[:,0] = ii
                    shift[:,1] = jj
                    shift[:,2] = kk
                    self.sc_basis_pos[count*self.p_num_atoms:
                            (count+1)*self.p_num_atoms,:] = self.p_basis_pos+shift
                    self.sc_basis_types.extend(self.p_basis_types)
                    count = count+1

        # reduced (crystal) coords
        self.sc_basis_pos[:,0] = self.sc_basis_pos[:,0]/self.sc_matrix[0,0]
        self.sc_basis_pos[:,1] = self.sc_basis_pos[:,1]/self.sc_matrix[1,1]
        self.sc_basis_pos[:,2] = self.sc_basis_pos[:,2]/self.sc_matrix[2,2]
                    
    # -----------------------------------------------------------------------------------------------

    def _minimum_image(self):

        """
        get relative position vectors (using minimum image) between basis atoms and supercell atoms
        and convert those to spherical coords too
        """

        # impose minimum image
        for atom in range(self.p_num_atoms):

            pos = self.sc_basis_pos[atom,:] # position of this atom in cartesian coords
            rel_pos = np.copy(self.sc_basis_pos) # relative positions
            rel_pos[:,0] = rel_pos[:,0]-pos[0]
            rel_pos[:,1] = rel_pos[:,1]-pos[1]
            rel_pos[:,2] = rel_pos[:,2]-pos[2]
            shift = np.zeros((self.sc_num_atoms,3))

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
            self.sc_rel_pos[atom,:,:] = rel_pos

            # convert relative position in crystal coords to cartesian
            self.sc_rel_cart[atom,:,:] = self._crystal_to_cart(self.sc_lat_vecs,self.sc_rel_pos[atom,:,:])

            # convert relative position from cartesian to spherical coords
            self.sc_rel_sph[atom,:,:] = self._cart_to_spherical(self.sc_rel_cart[atom,:,:])

# ===================================================================================================
# ---------------------------------------------------------------------------------------------------
# ===================================================================================================














