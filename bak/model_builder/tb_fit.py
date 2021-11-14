# system modules
import sys

# set path to custom modules
module_path = './modules'
sys.path.append(module_path)

# custom modules 
import m_lattice
import m_orbitals
import m_electron_hamiltonian

# --------------------------------------------------------------------------------------------------
# build a model hamiltonian

def model_init():

    """
    set up and return a model to do stuff with
    """

    # ----------------------------------------------------------------------------------------------
    # set up the crystal lattice and define atom positions 

    lattice = m_lattice.lattice()

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

    # ----------------------------------------------------------------------------------------------
    # setup supercell

    sc_scale = 3.929
    sc_vecs = [[2.000, 0.000, 0.000], # row vectors in angstroms
               [0.000, 2.000, 0.000],
               [0.000, 0.000, 2.496]]
    lattice.set_supercell(sc_vecs,sc_scale)

    # ----------------------------------------------------------------------------------------------
    # set up class to hold all the orbitals etc.

    orbitals = m_orbitals.orbitals(lattice)

    # now set the orbtials
    orbitals.set_orbital(site=0,n=3,orb='d_x2-y2')
    orbitals.set_orbital(site=0,n=3,orb='d_x2-y2')
    orbitals.set_orbital(site=1,n=2,orb='px')
    orbitals.set_orbital(site=1,n=2,orb='py')
    orbitals.set_orbital(site=2,n=2,orb='px')
    orbitals.set_orbital(site=2,n=2,orb='py')

    # define irreducible set of hoppings between orbitals and onsite energies 
    t0 = -2   # Cu-O 1st-nn hopping between Cu:3d_x2-y2 and O:2_p
    t1 = -0.5 # Cu-O 2nd-nn hopping
    t3 = -1   # O-O 1st-nn hopping
    t4 = -0.1 # O-O 2nd-nn hopping
    s0 = 1    # Cu onsite energy
    s1 = 0.2  # O onsite energy

    # now set the orbtials. set ALL of them using the irreducible set above
    orbtials.set_hopping()

    # ----------------------------------------------------------------------------------------------
    # setup hamiltonian

    eletr_hamlt = m_electron_hamiltonian.hamiltonian(lattice,orbitals)

    # ----------------------------------------------------------------------------------------------
    # return the model to be used for stuff

    return hamiltonian


# --------------------------------------------------------------------------------------------------
# if not imported as a module, do the stuff above

if __name__ == '__main__':

    model = model_init()









