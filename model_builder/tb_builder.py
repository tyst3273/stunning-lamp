# system modules

# custom modules 
import mod_lattice
import mod_hamiltonian

# --------------------------------------------------------------------------------------------------
# build a model hamiltonian

def model_init():

    """
    set up and return a model to do stuff with
    """
    
    # set up the crystal lattice and define atom positions 
    lattice = mod_lattice.lattice()

    # set primitive vectors
    prim_scale = 3.929
    prim_vecs = [[1.000, 0.000, 0.000], # row vectors in angstroms
                 [0.000, 1.000, 0.000],
                 [0.000, 0.000, 2.496]]
    lattice.set_prim_vecs(prim_vecs,prim_scale)

    # set basis atoms 
    basis_pos = [[0.0, 0.0, 0.0], # row vectors in reduced (crystal) coordinates
                 [0.5, 0.0, 0.5],
                 [0.0, 0.5, 0.0]]
    basis_labels = ['Cu','O','O'] # can be anything but will be used to identify inequivalent sites
    lattice.set_basis_positions(basis_pos,basis_labels)

    # setup supercell
    sc_scale = 3.929
    sc_vecs = [[2.000, 0.000, 0.000], # row vectors in angstroms
               [0.000, 2.000, 0.000],
               [0.000, 0.000, 2.496]]
    lattice.set_supercell(sc_vecs,sc_scale)


    # not as sure what to do here. need to set 'orbitals' on each site and define the coupling
    # between them. dunno how i want to do it.  


    # setup hamiltonian
    hamiltonian = mod_hamiltonian.hamiltonian(lattice)

    return hamiltonian

# --------------------------------------------------------------------------------------------------
# if not imported as a module, do the stuff above

if __name__ == '__main__':

    model = model_init()









