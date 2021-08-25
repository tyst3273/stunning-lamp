# system modules

# custom modules
import m_lattice
import m_atoms

def hamiltonian_init():

    # ----------------------------------------------------------------------------------------------
    # set up the lattice

    lattice = m_lattice.lattice()
    
    lat_vecs = [[0.0, 1.0, 1.0], # row vectors in angstrom
                [1.0, 0.0, 1.0],
                [1.0, 1.0, 0.0]]
    lat_scale = 1.98918 # lattice vector scaling factor
    lattice.set_lattice_vectors(lat_vecs,lat_scale)

    # ----------------------------------------------------------------------------------------------
    # set atom positions

    atoms = m_atoms.atoms()

    atom_pos = [[0.0, 0.0, 0.0]]


if __name__ == '__main__':

    hamlt = hamiltonian_init()




