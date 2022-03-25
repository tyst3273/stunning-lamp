# system modules
import numpy as np
#import spglib

# custom modules


# --------------------------------------------------------------------------------------------------

class crystal:

    def __init__(self,atom_pos,atom_types,lattice_vectors,cartesian=False):

        self.set_lattice_vectors(lattice_vectors)
        self.set_atom_positions(atom_pos,cartesian)
        self.set_atom_types(atom_types)
    
    # ----------------------------------------------------------------------------------------------

    def set_lattice_vectors(self,lattice_vectors):
        """
        lattice_vectors is a 3x3 numpy float array
        """
        try:
            lattice_vectors = np.array(lattice_vectors,dtype=float)
        except:
            exit('lattice_vectors should be list of list of floats, scale should be a float')
        if len(lattice_vectors.shape) != 2 or lattice_vectors.shape[0] != 3 or \
                lattice_vectors.shape[1] != 3:
            exit('lattice_vectors should be 3x3 (3D vectors)')
        self._lattice_vectors = lattice_vectors

    def set_atom_positions(self,atom_pos,cartesian):
        """
        atom_pos is an Nx3 numpy float array. cartesian is boolean
        """
        try:
            atom_pos = np.array(atom_pos,dtype=float)
        except:
            exit('atom_pos should be list of list of floats')
        if len(atom_pos.shape) != 2 or atom_pos.shape[1] != 3:
            exit('atom_pos should be Nx3 (3D vectors)')
        self._num_atoms = atom_pos.shape[0]
        if cartesian:
            self._cart_pos = atom_pos
            self._red_pos = self._cartesian_to_crystal_coords(self._lattice_vetors,self._cart_pos)
        else:
            self._red_pos = atom_pos
            self._cart_pos = self._crystal_to_cartesian_coords(self._lattice_vectors,self._red_pos)

    def set_atom_types(self,atom_types):
        """
        atom_types is an N element list of strings. it's convereted to numpy object array
        for quick manipulation
        """
        try:
            atom_types = np.array(atom_types,dtype=object).flatten()
        except:
            exit('atom_types should be a list of str')
        if atom_types.size != self._num_atoms:
            exit('number of atom_types doesnt match number of positions')
        self._unique_types, self._type_inds, self._type_counts = np.unique(atom_types,
                                    return_inverse=True,return_counts=True)
        self._num_types = self._unique_types.size
        self._atom_types = atom_types

    # ----------------------------------------------------------------------------------------------

    def _crystal_to_cartesian_coords(self,mat,vecs):
        return self._change_basis(mat,vecs)

    def _cartesian_to_crystal_coords(self,mat,vecs):
        mat = np.linalg.inv(mat)
        return self._change_basis(mat,vecs) 

    def _change_basis(self,mat,vecs):
        """
        lattice vectors is matrix of row vectors, but we need column vecs to change basis 
        correctly. so transpose.
        """
        mat = mat.T 
        pos = np.zeros(vecs.shape)
        for ii in range(vecs.shape[0]):
            pos[ii,:] = np.matmul(mat,vecs[ii,:])
        return pos

    # ----------------------------------------------------------------------------------------------

    def make_supercell(self,sc_matrix):
        """
        only "diagonal" supercell matrix allowed for now
        """
        try:
            self._sc_matrix = np.array(sc_matrix,dtype=int)
        except:
            exit('sc_matrix must be a list of 3 integers')

        self._num_formula_units = self._sc_matrix.prod()
        self._num_sc_atoms = self._num_formula_units*self._num_atoms

        self._sc_lattice_vectors = np.zeros((3,3))
        for ii in range(3):
            self._sc_lattice_vectors[ii,:] = self._lattice_vectors[ii,:]*self._sc_matrix[ii]

        self._sc_unitcell_inds = np.zeros(self._num_sc_atoms,dtype=int)
        self._sc_red_pos = np.zeros((self._num_sc_atoms,3))
        self._sc_types = np.zeros(self._num_sc_atoms,dtype=object)
        shift = 0
        cell = 0
        for ii in range(self._sc_matrix[0]):
            for jj in range(self._sc_matrix[1]):
                for kk in range(self._sc_matrix[2]):
                    self._sc_red_pos[shift:shift+self._num_atoms,:] = self._red_pos[:,:]
                    self._sc_red_pos[shift:shift+self._num_atoms,0] += ii
                    self._sc_red_pos[shift:shift+self._num_atoms,1] += jj
                    self._sc_red_pos[shift:shift+self._num_atoms,2] += kk
                    self._sc_types[shift:shift+self._num_atoms] = self._atom_types[:]
                    self._sc_unitcell_inds[shift:shift+self._num_atoms] = cell
                    cell += 1
                    shift += self._num_atoms

        self._sc_red_pos = self._sc_red_pos/self._sc_matrix
        self._sc_unique_types, self._sc_type_inds, self._sc_type_counts = np.unique( 
                        self._sc_types,return_inverse=True,return_counts=True)
        self._sc_cart_pos = self._crystal_to_cartesian_coords(self._sc_lattice_vectors,
                        self._sc_red_pos)

    # ----------------------------------------------------------------------------------------------

    def freeze_rotation_in_supercell(self,coord=[0,0,0],nn_dist=1,euler_angles=[12,0,0]):
        """
        coord is coordinate (in crystal units!) to center rotation around.
        pick all atoms within nn_dist (in angstroms) and rotate them
        according to Euler angles. Euler angles are in degrees. See goldstein sec 4.4
        for conventions
        """
        if not hasattr(self,'_sc_red_pos'):
            exit('no _sc_red_pos. did you use make_supercell?')
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
        rel, shift = self._get_minimum_image_arrays(coord,self._sc_red_pos)
        min_im_pos = rel+shift
        min_im_cart = self._crystal_to_cartesian_coords(self._sc_lattice_vectors,min_im_pos)
        min_im_dist = np.sqrt((min_im_cart**2).sum(axis=1))
        self._which_to_rotate = (min_im_dist <= nn_dist)
        nn_inds = np.argwhere(self._which_to_rotate).flatten()

        # apply the rotation
        for ii in range(nn_inds.size):
            ind = nn_inds[ii]
            min_im_cart[ind,:] = np.matmul(R_full,min_im_cart[ind,:])

        # now put pos back into crystal coords and unshift minimum image
        min_im_pos = self._cartesian_to_crystal_coords(self._sc_lattice_vectors,min_im_cart)
        min_im_pos = min_im_pos+coord-shift

        # replace sc positions with new ones
        self._sc_cart_pos = self._crystal_to_cartesian_coords(self._sc_lattice_vectors,min_im_pos)
        self._sc_red_pos = min_im_pos

    # ----------------------------------------------------------------------------------------------

    def _get_minimum_image_arrays(self,coord,pos):
        rel = pos-coord
        shift = (rel < -0.5).astype(int)-(rel > 0.5).astype(int)
        return rel, shift

    # ----------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------

class interface_io:

    """
    read and write from DFT calculators
    """

    def write_abivars(self,crystal,out_file='abivars',which='primitive'):
        if which == 'primitive':
            pos, lat_vecs, types, counts, inds = self._get_prim(crystal)
        elif which == 'supercell':
            pos, lat_vecs, types, counts, inds = self._get_sc(crystal)

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

    def write_poscar(self,crystal,out_file='POSCAR',which='primitive'):
        if which == 'primitive':
            pos, lat_vecs, types, counts, inds = self._get_prim(crystal)
        elif which == 'supercell':
            pos, lat_vecs, types, counts, inds = self._get_sc(crystal)

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
            for ii in range(crystal._num_types):
                _ += f' {types[ii]}'
            f_out.write(_)
            _ = '\n '
            for ii in range(crystal._num_types):
                _ += f' {counts[ii]}'
            f_out.write(_)
            f_out.write('\nDirect\n')
            for ii in range(num_atoms):
                f_out.write(f' {pos[ii,0]: 12.9f} {pos[ii,1]: 12.9f}' \
                            f' {pos[ii,2]: 12.9f} {types[inds[ii]]}\n')

    def write_xyz(self,crystal,out_file='pos.xyz',which='primitive'):
        if which == 'primitive':
            cart, lat_vecs, types, counts, inds = self._get_prim(crystal,cartesian=True)
        elif which == 'supercell':
            cart, lat_vecs, types, counts, inds = self._get_sc(crystal,cartesian=True)

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

    def write_cp2k(self,crystal,out_file='cp2k.pos',which='primitive'):
        if which == 'primitive':
            pos, lat_vecs, types, counts, inds = self._get_prim(crystal)
        elif which == 'supercell':
            pos, lat_vecs, types, counts, inds = self._get_sc(crystal)

        sort = np.argsort(inds).flatten()
        inds = inds[sort]
        pos = pos[sort,:]
        num_atoms = pos.shape[0]
        tab = ' '*4

        if hasattr(crystal,'_which_to_rotate'):
            rot = np.argwhere(crystal._which_to_rotate[sort]).flatten()
            _ = 'the following atoms have been \'rotated\':\n\t'
            for ii in range(rot.size):
                _ = _+f'{rot[ii]+1:g} '
            print(_)

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

    def _get_prim(self,crystal,cartesian=False):
        if cartesian:
            pos = crystal._cart_pos
        else:
            pos = crystal._red_pos
        lat_vecs = crystal._lattice_vectors
        types = crystal._unique_types
        counts = crystal._type_counts
        inds = crystal._type_inds
        return pos, lat_vecs, types, counts, inds

    def _get_sc(self,crystal,cartesian=False):
        if not hasattr(crystal,'_sc_red_pos'):
            exit('no _sc_red_pos. did you use make_supercell?')
        if cartesian:
            pos = crystal._sc_cart_pos
        else:
            pos = crystal._sc_red_pos
        lat_vecs = crystal._sc_lattice_vectors
        types = crystal._sc_unique_types
        counts = crystal._sc_type_counts
        inds = crystal._sc_type_inds
        return pos, lat_vecs, types, counts, inds

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





