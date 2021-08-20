# system modules
import numpy as np

# custom modules
from mod_utils import raise_error, print_message

# --------------------------------------------------------------------------------------------------
# (electronic) hamiltonian class

class hamiltonian:

    """
    hamiltonian takes lattice, atoms, bonds, etc. as arguments and sets up hamiltonian matrix at 
    whatever requested k-points. returns hamiltonian or diagonalizes it and returns eigenvals, 
    eigenvectors etc.
    """

    def __init__(self,lattice):

        """
        save the classes needed to build the hamiltonian
        """
        
        self.lattice = lattice

    # ----------------------------------------------------------------------------------------------

    def get_hamiltonian(self,k_point=[0,0,0]):

        """
        return the hamiltonian matrix at a given k-point
        """

        pass

    # ----------------------------------------------------------------------------------------------

    def get_eigenvals(self,k_point=[0,0,0]):

        """
        return the eigenvalues and eigenvectors at a given k-point
        """

        ham_k = self.get_hamitonian(k_point)
        eigvals, eigvecs = np.linalg.eigh(ham_k)
        
    # ----------------------------------------------------------------------------------------------




