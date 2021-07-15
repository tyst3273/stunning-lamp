import numpy as np
from mod_utils import print_stdout, raise_error


class lattice:

    # ------------------------------------------------------------------------------------------------
        
    def __init__(self):

        pass
        
    # ------------------------------------------------------------------------------------------------

    def set_atom_data(self,atom_masses):

        """
        set data for atoms. 
        """
        self.atom_masses = np.array(atom_masses) # atomic mass in amu

    # ------------------------------------------------------------------------------------------------

    def set_primitive_cell(self,prim_cell_file='prim_cell'):

        """
        parse file with atom data in it and set up primitive cell
        """

        # parse the file
        lat, nat, pos, types = self._parse_cell_file(prim_cell_file)

        # set the data
        self.prim_vecs = lat
        self.num_basis = nat
        self.basis_pos = pos
        self.basis_types = types

        # get the reciprocal lattice vectors
        self.r_prim_vecs, self.prim_vol = self._compute_reciprocal_lattice(self.prim_vecs)

        # get basis pos in cartesian coords
        self.basis_pos_cart = self._crystal_to_cart(self.prim_vecs,self.basis_pos)

        # check number of types mathces number of masses and number of atomic numbers
        self.num_types = np.unique(self.basis_types).shape[0]
        if self.num_types != self.atom_masses.shape[0]:
            message = 'number of basis types doesnt match number of masses'
            raise_error(message)

    # ------------------------------------------------------------------------------------------------

    def set_super_cell(self,super_cell_file='super_cell'):

        """
        parse file with atom data in it and set up super cell
        """

        # parse the file
        lat, nat, pos, types = self._parse_cell_file(super_cell_file)

        # set the data
        self.sc_vecs = lat
        self.num_sc = nat
        self.sc_pos = pos
        self.sc_types = types

        # get the reciprocal lattice vectors
        self.r_sc_vecs, self.sc_vol = self._compute_reciprocal_lattice(self.sc_vecs)

        # get sc pos in cartesian coords
        self.sc_pos_cart = self._crystal_to_cart(self.sc_vecs,self.sc_pos)

        # check number of types mathces number of masses and number of atomic numbers
        num_types = np.unique(self.sc_types).shape[0]
        if num_types != self.atom_masses.shape[0]:
            message = 'number of sc types doesnt match number of masses'
            raise_error(message)

    # -----------------------------------------------------------------------------------------------

    def _parse_cell_file(self,cell_file):

        """
        read a cell file
        """

        # read file
        with open(cell_file,'r') as f_cell:
            lines = f_cell.readlines()

        # get scale of lattice vectors
        try:
            lat_scale = float(lines[0].strip())
        except:
            message = f'the lattice scale factor in  \'{cell_file}\' seems wrong'
            raise_error(message)

        # get lattice vectors
        lat = np.zeros((3,3))
        try:
            for ii in range(3):
                vec = lines[ii+1].strip().split()
                lat[ii,:] = [float(x) for x in vec]
        except:
            message = f'the lattice vectors in \'{cell_file}\' seem wrong'
            raise_error(message)

        # get number of atoms to read
        try:
            nat = int(lines[4].strip())
        except:
            message = f'number of atoms line in \'{cell_file}\' should be an integer'
            raise_error(message)

        # get positions and types
        pos = np.zeros((nat,3))
        types = np.zeros(nat,dtype=int)
        try:
            for ii in range(nat):
                line = lines[ii+5].strip().split()
                types[ii] = line[0]
                pos[ii,:] = [float(x) for x in line[1:]]
        except:
            message = f'positions or type data in \'{cell_file}\' seem wrong'
            raise_error(message)

        # return it all
        return lat*lat_scale, nat, pos, types

    # -------------------------------------------------------------------------------------------------

    def find_neighbors(self):

        """
        get relative position vectors (using minimum image) between atoms in supercell
        and convert those to cartesian and spherical coords
        """
        # relative pos vectors. same vector, different coordinate basis
        self.sc_rel = np.zeros((self.num_sc,self.num_sc,3)) # crystal coords
        self.sc_rel_cart = np.zeros((self.num_sc,self.num_sc,3)) # cartesian coords
        self.sc_rel_sph = np.zeros((self.num_sc,self.num_sc,3)) # spherical coords

        for atom in range(self.num_sc):

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
            self.sc_rel[atom,:,:] = rel_pos

            # convert relative position in crystal coords to cartesian
            self.sc_rel_cart[atom,:,:] = self._crystal_to_cart(self.sc_vecs,self.sc_rel[atom,:,:])

            # convert relative position from cartesian to spherical coords
            self.sc_rel_sph[atom,:,:] = self._cart_to_spherical(self.sc_rel_cart[atom,:,:])

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
        cart_pos[ind,2] = 1
        sph_coords[:,1] = np.arccos(cart_pos[:,2]/r)  # theta

        # atan is defined for arg -> inf
        sph_coords[:,2] = np.arctan2(cart_pos[:,1],cart_pos[:,0])   # phi

        return sph_coords

    # ------------------------------------------------------------------------------------------------

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

    def _crystal_to_cart(self,lat_vecs,crystal_coords):

        """
        wrapper to just matrix multiply each atoms vector by the fiven transformation matrix.
        """

        nat = crystal_coords.shape[0]
        cart_coords = np.zeros((nat,3))
        for atom in range(nat):
            cart_coords[atom,:] = np.matmul(lat_vecs,crystal_coords[atom,:])

        return cart_coords

    # -----------------------------------------------------------------------------------------------
    # ---------------------------------- write info to files ----------------------------------------
    # -----------------------------------------------------------------------------------------------

    def write_lattice(self):

        """
        print all of the stuff setup in this class to screen for diagnostics
        """

        fout = open('LATTICE.OUT','w')

        # print primitive cell
        message = ' ** PRIMITIVE CELL **\n real space primitive cell (Angstrom):\n'
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
        fout.write(message)

        # print supercell
        message = ' ** SUPER CELL **\n real space supercell (Angstrom):\n'
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
        fout.write(message)

        fout.close()

    # -----------------------------------------------------------------------------------------------

    def write_neighbors(self):

        """
        print neighbor vectors to screen for diagnostics
        """
        
        fout = open('NEIGHBORS.OUT','w')

        # print relative pos
        message = f' ** NEIGHBORS **\n there are {self.num_basis} atoms in the primitive cell\n'
        for ii in range(self.num_sc):
            message = message+f'  sc atom index: {ii},  type: {self.sc_types[ii]},'\
                    f' pos (Angstrom): {self.sc_pos_cart[ii,0]: 2.3f}'\
                    f' {self.sc_pos_cart[ii,1]: 2.3f} {self.sc_pos_cart[ii,2]: 2.3f}\n'
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
        fout.write(message)

        fout.close()

    # -----------------------------------------------------------------------------------------------

# ===================================================================================================
# ---------------------------------------------------------------------------------------------------
# ===================================================================================================
   



