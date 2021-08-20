# system modules
import numpy as np

# custom modules
from m_utils import raise_error, print_message

# --------------------------------------------------------------------------------------------------
# (electronic) hamiltonian class

class hamiltonian:

    """
    hamiltonian takes lattice, atoms, bonds, etc. as arguments and sets up hamiltonian matrix at 
    whatever requested k-points. returns hamiltonian or diagonalizes it and returns eigenvals, 
    eigenvectors etc.
    """

    def __init__(self,lattice,orbitals):

        """
        save the classes needed to build the hamiltonian
        """
        
        self.lattice = lattice
        self.orbitals = orbitals

    # ----------------------------------------------------------------------------------------------

    def get_hamiltonian_k(self,k_point=[0,0,0]):

        """
        return the hamiltonian matrix at a given k-point
        """

        # this should do the space FT and return the hamiltonian at a given k-point
        # do something like this
        #ham_k = self._do_space_ft(k_point)
        #return ham_k

        pass

    # ----------------------------------------------------------------------------------------------

    def get_eigenvals_k(self,k_point=[0,0,0]):

        """
        return the eigenvalues and eigenvectors at a given k-point
        """

        # do something like this
        #ham_k = self.get_hamitonian_k(k_point)
        #eigvals, eigvecs = np.linalg.eigh(ham_k)
        
        pass

    # ----------------------------------------------------------------------------------------------




