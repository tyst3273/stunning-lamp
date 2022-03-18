# system modules
import numpy as np
#import spglib

# custom modules


# --------------------------------------------------------------------------------------------------

class crystal:

    def __init__(self,basis_pos,basis_types,lattice_vectors,magnetic_moments=None):

        self.lattice_vectors = lattice_vectors
        self.basis_pos = basis_pos
        self.basis_types = basis_types
        self.magnetic_moments = magnetic_moments
        self.io = interface_io()
    
    # ----------------------------------------------------------------------------------------------

    @property
    def lattice_vectors(self):
        return self._lattice_vectors
    @lattice_vectors.setter
    def lattice_vectors(self,value):
        try:
            value = np.array(value,dtype=float)
        except:
            exit('lattice_vectors should be list of list of floats')
        if len(value.shape) != 2 or value.shape[0] != 3 or value.shape[1] != 3:
            exit('lattice_vectors should be 3x3 (3D vectors)')
        #self._num_basis_atoms = value.shape[0]
        self._lattice_vectors = value

    @property
    def basis_pos(self):
        return self._basis_pos
    @basis_pos.setter
    def basis_pos(self,value):
        try:
            value = np.array(value,dtype=float)
        except:
            exit('basis_pos should be list of list of floats')
        if len(value.shape) != 2 or value.shape[1] != 3:
            exit('basis_pos should be Nx3 (3D vectors)')
        self.num_basis_atoms = value.shape[0]
        self._basis_pos = value
        self.basis_cart = self._crystal_to_cartesian_coords(self.lattice_vectors,self.basis_pos)

    @property
    def basis_types(self):
        return self._basis_types
    @basis_types.setter
    def basis_types(self,value):
        try:
            value = np.array(value,dtype=object).flatten()
        except:
            exit('basis_types should be a list of str')
        if value.size != self.num_basis_atoms:
            exit('number of basis_types doesnt match number of positions')
        self.unique_types, self.type_inds, self.type_counts = np.unique(value,
                                    return_inverse=True,return_counts=True)
        self.num_basis_types = self.unique_types.size
        self._basis_types = value

    @property
    def magnetic_moments(self):
        return self._magnetic_moments
    @magnetic_moments.setter
    def magnetic_moments(self,value):
        self._magnetic_moments = value

    # ----------------------------------------------------------------------------------------------

    def _crystal_to_cartesian_coords(self,mat,vecs):
        """
        note, need mat in col vectors not row, so transpose
        """
        mat = mat.T # need col vectors, not row vectors
        cart = np.zeros((vecs.shape[0],3))
        for ii in range(vecs.shape[0]):
            cart[ii,:] = np.matmul(mat,vecs[ii,:])
        return cart

    def _cartesian_to_crystal_coords(self,mat,vecs):
        """
        note, need mat in col vectors not row, so transpose
        to go from cart. to crystal, use inverse of lattice vectors
        """
        mat = np.linalg.inv(mat).T # need col vectors, not row vectors
        pos = np.zeros((vecs.shape[0],3))
        for ii in range(vecs.shape[0]):
            pos[ii,:] = np.matmul(mat,vecs[ii,:])
        return pos

    # ----------------------------------------------------------------------------------------------

    def make_supercell(self,sc_matrix):
        """
        only "diagonal" supercell matrix allowed for now
        """
        try:
            self.sc_matrix = np.array(sc_matrix,dtype=int)
        except:
            exit('sc_matrix must be a list of 3 integers')

        self.num_formula_units = self.sc_matrix.prod()
        self.num_sc_atoms = self.num_formula_units*self.num_basis_atoms

        self.sc_lattice_vectors = np.zeros((3,3))
        for ii in range(3):
            self.sc_lattice_vectors[ii,:] = self.lattice_vectors[ii,:]*self.sc_matrix[ii]

        self.sc_unitcell_inds = np.zeros(self.num_sc_atoms,dtype=int)
        self.sc_pos = np.zeros((self.num_sc_atoms,3))
        self.sc_types = np.zeros(self.num_sc_atoms,dtype=object)
        shift = 0
        cell = 0
        for ii in range(self.sc_matrix[0]):
            for jj in range(self.sc_matrix[1]):
                for kk in range(self.sc_matrix[2]):
                    self.sc_pos[shift:shift+self.num_basis_atoms,:] = self.basis_pos[:,:]
                    self.sc_pos[shift:shift+self.num_basis_atoms,0] += ii
                    self.sc_pos[shift:shift+self.num_basis_atoms,1] += jj
                    self.sc_pos[shift:shift+self.num_basis_atoms,2] += kk
                    self.sc_types[shift:shift+self.num_basis_atoms] = self.basis_types[:]
                    self.sc_unitcell_inds[shift:shift+self.num_basis_atoms] = cell
                    cell += 1
                    shift += self.num_basis_atoms

        self.sc_pos = self.sc_pos/self.sc_matrix
        self.sc_unique_types, self.sc_type_inds, self.sc_type_counts = np.unique( 
                        self.sc_types,return_inverse=True,return_counts=True)
        self.sc_cart = self._crystal_to_cartesian_coords(self.sc_lattice_vectors,self.sc_pos)

    # ----------------------------------------------------------------------------------------------

    def freeze_rotation_in_supercell(self,coord=[0,0,0],nn_dist=1,euler_angles=[12,0,0]):
        """
        coord is coordinate (in crystal units!) to center rotation around.
        pick all atoms within nn_dist (in angstroms) and rotate them
        according to Euler angles. Euler angles are in degrees. See goldstein sec 4.4
        for conventions
        """
        if not hasattr(self,'sc_pos'):
            exit('no sc_pos. did you use make_supercell?')
        try:
            euler_angles = np.array(euler_angles,dtype=float).flatten()
        except:
            exit('euler_angles seem wrong in freeze_rotation_in_supercell')
        if euler_angles.size != 3:
            exit('euler_angles seem wrong in freeze_rotation_in_supercell')

        # setup rotation matricies
        phi = euler_angles[0]*np.pi/180
        R_phi = np.array([[  np.cos(phi), -np.sin(phi),  0.0],
                          [  np.sin(phi),  np.cos(phi),  0.0],
                          [         0.0,           0.0,  1.0]])
        theta = euler_angles[1]*np.pi/180
        R_theta = np.array([[  1.0,            0.0,            0.0],
                            [  0.0,  np.cos(theta), -np.sin(theta)],
                            [  0.0,  np.sin(theta),  np.cos(theta)]])
        psi = euler_angles[2]*np.pi/180
        R_psi = np.array([[ np.cos(psi), -np.sin(psi),  0.0],
                          [ np.sin(psi),  np.cos(psi),  0.0],
                          [         0.0,          0.0,  1.0]])
        R_full = np.matmul(R_psi,np.matmul(R_theta,R_phi))

        # find nn
        rel, shift = self._get_minimum_image_arrays(coord,self.sc_pos)
        min_im_pos = rel+shift
        min_im_cart = self._crystal_to_cartesian_coords(self.sc_lattice_vectors,min_im_pos)
        min_im_dist = np.sqrt((min_im_cart**2).sum(axis=1))
        nn_inds = np.argwhere(min_im_dist <= nn_dist).flatten()

        # apply the rotation
        for ii in range(nn_inds.size):
            ind = nn_inds[ii]
            min_im_cart[ind,:] = np.matmul(R_full,min_im_cart[ind,:])

        # now put pos back into crystal coords and unshift minimum image
        min_im_pos = self._cartesian_to_crystal_coords(self.sc_lattice_vectors,min_im_cart)
        min_im_pos = min_im_pos+coord-shift

        # replace sc positions with new ones
        self.sc_cart = self._crystal_to_cartesian_coords(self.sc_lattice_vectors,min_im_pos)
        self.sc_pos = min_im_pos

    # ----------------------------------------------------------------------------------------------

    def _get_minimum_image_arrays(self,coord,pos):
        rel = pos-coord
        shift = (rel < -0.5).astype(int)-(rel > 0.5).astype(int)
        return rel, shift

    # ----------------------------------------------------------------------------------------------

    def write_poscar(self,out_file='POSCAR',which='primitive'):
        self.io.write_poscar(self,out_file,which)

    def write_xyz(self,out_file='pos.xyz',which='primitive'):
        self.io.write_xyz(self,out_file,which)

    def write_cp2k(self,out_file='cp2k.pos',which='primitive'):
        self.io.write_cp2k(self,out_file,which)

    def write_abivars(self,out_file='abivars',which='primitive'):
        self.io.write_abivars(self,out_file,which)

    # ----------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------

class interface_io:

    """
    read and write from DFT calculators
    """

    def write_abivars(self,crystal,out_file,which):
        if which == 'primitive':
            pos = crystal.basis_pos
            lat_vecs = crystal.lattice_vectors
            types = crystal.unique_types
            counts = crystal.type_counts
            inds = crystal.type_inds
        elif which == 'supercell':
            if not hasattr(crystal,'sc_pos'):
                exit('no sc_pos. did you use make_supercell?')
            pos = crystal.sc_pos
            lat_vecs = crystal.sc_lattice_vectors
            types = crystal.sc_unique_types
            counts = crystal.sc_type_counts
            inds = crystal.sc_type_inds

        sort = np.argsort(inds).flatten()
        inds = inds[sort]
        pos = pos[sort,:]
        num_atoms = pos.shape[0]
        with open(out_file,'w') as f_out:
            f_out.write('# Created by my custom scripts\n')
            f_out.write('\nacell  3*1.00 Angstrom\n')
            f_out.write(f'\nrprim {lat_vecs[0,0]: 12.9f} {lat_vecs[0,1]: 12.9f}'\
                             f' {lat_vecs[0,2]: 12.9f}\n')
            f_out.write(f'      {lat_vecs[1,0]: 12.9f} {lat_vecs[1,1]: 12.9f}'\
                             f' {lat_vecs[1,2]: 12.9f}\n')
            f_out.write(f'      {lat_vecs[2,0]: 12.9f} {lat_vecs[2,1]: 12.9f}'\
                             f' {lat_vecs[2,2]: 12.9f}\n')
            f_out.write(f'\nnatom  {num_atoms}\n')
            f_out.write('\nxred_symbols\n')
            for ii in range(num_atoms):
                f_out.write(f' {pos[ii,0]: 12.9f} {pos[ii,1]: 12.9f}' \
                            f' {pos[ii,2]: 12.9f} {types[inds[ii]]}\n')

    def write_poscar(self,crystal,out_file,which):
        if which == 'primitive':
            pos = crystal.basis_pos
            lat_vecs = crystal.lattice_vectors
            types = crystal.unique_types
            counts = crystal.type_counts
            inds = crystal.type_inds
        elif which == 'supercell':
            if not hasattr(crystal,'sc_pos'):
                exit('no sc_pos. did you use make_supercell?')
            pos = crystal.sc_pos
            lat_vecs = crystal.sc_lattice_vectors
            types = crystal.sc_unique_types
            counts = crystal.sc_type_counts
            inds = crystal.sc_type_inds

        sort = np.argsort(inds).flatten()
        inds = inds[sort]
        pos = pos[sort,:]
        num_atoms = pos.shape[0]
        with open(out_file,'w') as f_out:
            f_out.write('Created by my custom scripts\n')
            f_out.write(' 1.00\n')
            for ii in range(3):
                f_out.write(f' {lat_vecs[ii,0]: 12.9f} {lat_vecs[ii,1]: 12.9f}'\
                             f' {lat_vecs[ii,2]: 12.9f}\n')            
            _ = ''
            for ii in range(crystal.num_basis_types):
                _ += f' {types[ii]}'
            f_out.write(_)
            _ = '\n '
            for ii in range(crystal.num_basis_types):
                _ += f' {counts[ii]}'
            f_out.write(_)
            f_out.write('\nDirect\n')
            for ii in range(num_atoms):
                f_out.write(f' {pos[ii,0]: 12.9f} {pos[ii,1]: 12.9f}' \
                            f' {pos[ii,2]: 12.9f} {types[inds[ii]]}\n')

    def write_xyz(self,crystal,out_file,which):
        if which == 'primitive':
            cart = crystal.basis_cart
            lat_vecs = crystal.lattice_vectors
            types = crystal.unique_types
            counts = crystal.type_counts
            inds = crystal.type_inds
        elif which == 'supercell':
            if not hasattr(crystal,'sc_pos'):
                exit('no sc_pos. did you use make_supercell?')
            cart = crystal.sc_cart
            lat_vecs = crystal.sc_lattice_vectors
            types = crystal.sc_unique_types
            counts = crystal.sc_type_counts
            inds = crystal.sc_type_inds

        sort = np.argsort(inds).flatten()
        inds = inds[sort]
        cart = cart[sort,:]
        num_atoms = cart.shape[0]
        with open(out_file,'w') as f_out:
            f_out.write(f' {num_atoms}\n')
            f_out.write('Created by my custom scripts\n')
            for ii in range(num_atoms):
                f_out.write(f' {types[inds[ii]]} {cart[ii,0]: 12.9f} {cart[ii,1]: 12.9f}' \
                            f' {cart[ii,2]: 12.9f}\n')

    def write_cp2k(self,crystal,out_file,which):
        if which == 'primitive':
            pos = crystal.basis_pos
            lat_vecs = crystal.lattice_vectors
            types = crystal.unique_types
            counts = crystal.type_counts
            inds = crystal.type_inds
        elif which == 'supercell':
            if not hasattr(crystal,'sc_pos'):
                exit('no sc_pos. did you use make_supercell?')
            pos = crystal.sc_pos
            lat_vecs = crystal.sc_lattice_vectors
            types = crystal.sc_unique_types
            counts = crystal.sc_type_counts
            inds = crystal.sc_type_inds

        sort = np.argsort(inds).flatten()
        inds = inds[sort]
        pos = pos[sort,:]
        num_atoms = pos.shape[0]
        tab = ' '*4
        with open(out_file,'w') as f_out:
            f_out.write(2*tab+f'&CELL\n')
            f_out.write(3*tab+f'A {lat_vecs[0,0]: 12.9f} {lat_vecs[0,1]: 12.9f}'\
                                f' {lat_vecs[0,2]: 12.9f}\n')
            f_out.write(3*tab+f'B {lat_vecs[1,0]: 12.9f} {lat_vecs[1,1]: 12.9f}'\
                                f' {lat_vecs[1,2]: 12.9f}\n')
            f_out.write(3*tab+f'C {lat_vecs[2,0]: 12.9f} {lat_vecs[2,1]: 12.9f}'\
                                f' {lat_vecs[2,2]: 12.9f}\n')
            f_out.write(3*tab+f'PERIODIC XYZ\n')
            f_out.write(2*tab+f'&END CELL\n')
            f_out.write(2*tab+f'&COORD\n')
            f_out.write(3*tab+f'SCALED\n')
            for ii in range(num_atoms):
                f_out.write(3*tab+f'{types[inds[ii]]:<2} {pos[ii,0]: 12.9f}'\
                                f' {pos[ii,1]: 12.9f} {pos[ii,2]: 12.9f}\n')
            f_out.write(2*tab+f'&END COORD')

    # ----------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------

class wave_vectors:

    def __init__(self,crystal):
        
        """
        create full and reduced grids of wave-vectors. use phonopy
        generate a kpt grid from the crystal symmetry (using) spglib
        have methods to write kpt grid files for different codes
        """
        pass

# --------------------------------------------------------------------------------------------------





