# system modules

# custom modules
from mod_model import model


# -------------------------------------------------------------------------------------------------
# return a model hamiltonian

def model_init():
    
    """
    sets up and returns a 'model' hamiltonian
    """

    tb_model = model()

    # set lattice vectors
    prim_scale = 3.929
    prim_vecs = [[1.000, 0.000, 0.000], # row vectors in angstroms
                 [0.000, 1.000, 0.000],
                 [0.000, 0.000, 2.496]]
    tb_model.set_lattice(prim_vecs,prim_scale)

    return tb_model


# -------------------------------------------------------------------------------------------------
# if not being called from another file, do the thing

if __name__ == '__main__':

    cuprate = model_init()








