import numpy as np
from mod_utils import print_stdout, raise_error
from mod_constants import hbar_sq_over_2m


# ====================================================================================================
# ----------------------------------------------------------------------------------------------------
# ====================================================================================================

class pw:

    # ------------------------------------------------------------------------------------------------

    def __init__(self):
        pass

    # ------------------------------------------------------------------------------------------------

    def set_plane_waves(self,crystal,pw_cut_off):

        # cut off should be a number
        try:
            pw_cut_off = float(pw_cut_off)
        except:
            message = 'pw_cutoff should a be a floating point number that is larger than 0'
            raise_error(message)

        # cut off should be larger than 0
        if pw_cut_off <= 0:
            message = 'pw_cutoff should a be a floating point number that is larger than 0'
            raise_error(message)

        self.pw_cut_off = pw_cut_off
       
        # find shortest recip. lat vec
        b0 = self._norm(crystal.r_lat_vecs[0,:])
        for ii in range(1,3):
            bi = self._norm(crystal.r_lat_vecs[ii,:])
            if bi <= b0:
                b0 = bi

        # find largest n to check
        n_max = int(1.5*np.sqrt(pw_cut_off/hbar_sq_over_2m)/b0)

        # create list of G to check
        num_test = (2*n_max+1)**3 # i.e. [-n_max,n_max]*[-n_max,n_max]*[-n_max,n_max]
        test_G = np.zeros((num_test,3))
        test_E = np.zeros(num_test)
    
        # loop over G vectors and find which ones are below cutoff
        ind = 0
        for nx in range(-n_max,n_max+1): # [-n_max,n_max]
            for ny in range(-n_max,n_max+1): # [-n_max,n_max]
                for nz in range(-n_max,n_max+1): # [-n_max,n_max]
                    test_G[ind,:] = [nx,ny,nz]
                    test_E[ind] = self._norm(nx*crystal.r_lat_vecs[0,:]+ny*crystal.r_lat_vecs[1,:]+
                                             nz*crystal.r_lat_vecs[2,:])**2 # G**2
                    ind = ind+1
        
        # find which are below cutoff
        test_E = test_E*hbar_sq_over_2m # KE = hbar**2/2/m * G**2
        test_E = np.argwhere(test_E <= pw_cut_off).flatten() # find inds below cutoff
        self.G_vecs = test_G[test_E] # PWs in basis set
        self.num_G = self.G_vecs.shape[0] # number of PWs in basis set

        # print number of plane waves
        message = f'there are {self.num_G} plane waves in the basis set'
        print_stdout(message,msg_type='PLANE WAVES')

    # ------------------------------------------------------------------------------------------------

    def _norm(self,vec):

        return np.sqrt(vec[0]**2+vec[1]**2+vec[2]**2)

    # ------------------------------------------------------------------------------------------------



