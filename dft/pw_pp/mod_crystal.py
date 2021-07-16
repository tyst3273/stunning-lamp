import numpy as np
from mod_utils import raise_error, print_stdout


# ====================================================================================================
# ----------------------------------------------------------------------------------------------------
# ====================================================================================================

class crystal:

    # -------------------------------------------------------------------------------------------------

    def __init__(self):
        pass

    # -------------------------------------------------------------------------------------------------

    def set_lattice(self,lat_vecs,lat_scale):

        # check that lattice vector matrix is correct dims
        self.lat_vecs = np.array(lat_vecs)
        if self.lat_vecs.shape[0] !=3 or self.lat_vecs.shape[1] !=3:
            message = 'lattice vectors should be a 3x3 matrix'
            raise_error(message)
        self.lat_vecs = self.lat_vecs*lat_scale

        # compute reciprocal lattice
        self._compute_reciprocal_lattice()

        # print lattice info
        message = 'real space lattice (Angstrom):\n'
        for ii in range(3):
            message = message+f'  {self.lat_vecs[ii,0]: 2.4f} {self.lat_vecs[ii,1]: 2.4f}' \
                    f' {self.lat_vecs[ii,2]: 2.4f}\n'
        message = message+' reciprocal lattice (1/Angstrom):\n'
        for ii in range(3):
            message = message+f'  {self.r_lat_vecs[ii,0]: 2.4f} {self.r_lat_vecs[ii,1]: 2.4f}' \
                    f' {self.r_lat_vecs[ii,2]: 2.4f}\n'
        print_stdout(message,msg_type='LATTICE')

    # -------------------------------------------------------------------------------------------------

    def _compute_reciprocal_lattice(self):

        self.r_lat_vecs = np.zeros((3,3))
        self.cell_vol = self.lat_vecs[0,:].dot(np.cross(self.lat_vecs[1,:],self.lat_vecs[2,:]))
        self.r_lat_vecs[0,:] = 2*np.pi*np.cross(self.lat_vecs[1,:],self.lat_vecs[2,:])/self.cell_vol
        self.r_lat_vecs[1,:] = 2*np.pi*np.cross(self.lat_vecs[2,:],self.lat_vecs[0,:])/self.cell_vol
        self.r_lat_vecs[2,:] = 2*np.pi*np.cross(self.lat_vecs[0,:],self.lat_vecs[1,:])/self.cell_vol

    # -------------------------------------------------------------------------------------------------

    def set_atoms(self,atom_pos,atom_types,atomic_numbers):

        # check atomic positions make sense
        self.atom_pos = np.array(atom_pos)
        if self.atom_pos.shape[1] != 3:
            message = 'atom positions should be in 3D'
            raise_error(message)

        # check that there is 1 type per atom
        self.num_atoms = self.atom_pos.shape[0]
        if self.num_atoms != len(atom_types):
            message = 'there shoulb be 1 type per atom'
            raise_error(message)
        self.atom_types = np.array(atom_types)

        # check there is 1 atomic number per type
        self.num_types = np.unique(self.atom_types).shape[0]
        self.atomic_numbers = np.array(atomic_numbers)
        if self.num_types != self.atomic_numbers.shape[0]:
            message = 'there should be one atomic number per atom type'
            raise_error(message)

        # convert crystal to cartesian coords
        self._transform_to_cartesian_coords()

        # print coords
        message = 'index, type, positions in crystal coords:\n'
        for atom in range(self.num_atoms):
            message = message+f'   {atom} {self.atom_types[atom]}   {self.atom_pos[atom,0]: 2.4f}' \
                              f' {self.atom_pos[atom,1]: 2.4f} {self.atom_pos[atom,2]: 2.4f}\n'
        message = message+' index, type, positions in cartesian coords (Angstrom):\n'
        for atom in range(self.num_atoms):
            message = message+f'   {atom} {self.atom_types[atom]}   {self.cart_pos[atom,0]: 2.4f}' \
                              f' {self.cart_pos[atom,1]: 2.4f} {self.cart_pos[atom,2]: 2.4f}\n'
        print_stdout(message,msg_type='COORDINATES')

    # -------------------------------------------------------------------------------------------------

    def _transform_to_cartesian_coords(self):

        self.cart_pos = np.zeros((self.num_atoms,3))
        for atom in range(self.num_atoms):
            self.cart_pos[atom,:] = np.matmul(self.lat_vecs,self.atom_pos[atom,:])

    # -------------------------------------------------------------------------------------------------


# =====================================================================================================
# -----------------------------------------------------------------------------------------------------
# =====================================================================================================



