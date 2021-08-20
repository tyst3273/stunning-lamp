# system modules

# custom modules
from m_utils import raise_error

# --------------------------------------------------------------------------------------------------
# orbitals class. should hold all the orbitals and hopping/onsite params

class orbitals:

    """
    hold the orbital definitions and hopping/onsite params
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self,lattice):

        """
        copy classes and define set of allowed orbitals
        """

        # copy the class(es) over
        self.lattice = lattice

        # a list of the allowed orbital types
        self.allowed_orbitals = ['s',                                   # the s orbital
                                 'p_x','p_y','p_z',                     # the p orbitals
                                 'd_z2','d_xz','d_yz','d_xy','d_x2-y2'] # the d orbitals

        # a dictionary that holds which orbital is on which site
        self.orbs = {} 

    # ----------------------------------------------------------------------------------------------

    def set_orbital(self,site,n,orb):

        # check that the requested orbital is allowed
        if orb not in allowed_orbitals:
            message = f'the orbital type {orb} is not allowed'
            raise_error(message)

        # check the the requested orbital doesnt already exist
        handle = f'site_{site}'
        if handle in self.orbs.keys():
            print(handle)
        else:
            self.orbs[handle]['site'] = site
            self.orbs[handle]['n'] = n
            self.orbs[handle]['orb'] = orb

    # ----------------------------------------------------------------------------------------------

    def set_hopping(self):

        # should set the hopping param between the 2 orbitals on 2 sites. they can be the same 
        # site as long as they are different {l,m}, e.g. like molecular hopping. should only allow
        # hopping between sites and orbitals that exist and should not allow hopping between the 
        # same pairs

        pass

    # ----------------------------------------------------------------------------------------------

