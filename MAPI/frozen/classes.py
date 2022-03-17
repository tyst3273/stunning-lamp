# system modules
import numpy as np
#import spglib

# custom modules


# --------------------------------------------------------------------------------------------------

class crystal:
    
    # ----------------------------------------------------------------------------------------------


    def __init__(self,basis_pos,basis_types,lattice_vectors,magnetic_moments=None):

        self.lattice_vectors = lattice_vectors
        self.basis_pos = basis_pos
        self.basis_types = basis_types
        self.magnetic_moments = magnetic_moments
    
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

    # ----------------------------------------------------------------------------------------------

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

    # ----------------------------------------------------------------------------------------------

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
        self._basis_types = value

    # ----------------------------------------------------------------------------------------------

    @property
    def magnetic_moments(self):
        return self._magnetic_moments
    @magnetic_moments.setter
    def magnetic_moments(self,value):
        self._magnetic_moments = value

    # ----------------------------------------------------------------------------------------------

    def _crystal_to_cartesian_coords(self,mat,vecs):
        """
        note, since lattice vectors are row but pos are 
        """
        cart = np.zeros((vecs.shape[0],3))
        for ii in range(vecs.shape[0]):
            cart[ii,:] = np.matmul(vecs[ii,:],mat)
        return cart

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

        self.sc_pos = np.zeros((self.num_sc_atoms,3))
        self.sc_types = np.zeros(self.num_sc_atoms,dtype=object)
        shift = 0
        for ii in range(self.sc_matrix[0]):
            for jj in range(self.sc_matrix[1]):
                for kk in range(self.sc_matrix[2]):
                    self.sc_pos[shift:shift+self.num_basis_atoms,:] = self.basis_pos[:,:]
                    self.sc_pos[shift:shift+self.num_basis_atoms,0] += ii
                    self.sc_pos[shift:shift+self.num_basis_atoms,1] += jj
                    self.sc_pos[shift:shift+self.num_basis_atoms,2] += kk
                    self.sc_types[shift:shift+self.num_basis_atoms] = self.basis_types[:]
                    shift += self.num_basis_atoms
        self.sc_pos = self.sc_pos/self.sc_matrix

        self.unique_sc_types, self.sc_type_inds, self.sc_type_counts = np.unique( 
                        self.sc_types,return_inverse=True,return_counts=True)
        self.sc_cart = self._crystal_to_cartesian_coords(self.sc_lattice_vectors,self.sc_pos)

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





