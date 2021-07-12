import numpy as np
from mod_utils import print_stdout, raise_error
from mod_constants import hbar_sq_over_2m


# ====================================================================================================
# ----------------------------------------------------------------------------------------------------
# ====================================================================================================

class kpoints:

    # ------------------------------------------------------------------------------------------------

    def __init__(self):
        
        self.allowed_reductions = ['full']

    # ------------------------------------------------------------------------------------------------

    def set_grid(self,kpt_grid,grid_reduction='full'):

        # check that kpoint grid makes sense
        try:
            self.kpt_grid = np.array(kpt_grid,dtype=int)
        except:
            message = 'the kpoint grid should be a list of 3 integers'
            raise_error(message)
        if self.kpt_grid.shape[0] != 3:
            message = 'the kpoint grid should be a list of 3 integers'
            raise_error(message)

        # check that requested method makes sense
        if grid_reduction not in self.allowed_reductions:
            message = f'kpoint grid reduction method \'{grid_reduction}\' not one of the' \
                    ' allowed methods.'
            raise_error(message)

        # grid reduction method
        self.grid_reduction = grid_reduction

        # set up grid based on requested method
        if self.grid_reduction == 'full':
            self._set_full_grid()

        # print the kpt grid
        message = f'there are {self.num_k} kpoints: (only printing the 1st 50)\n'
        message = message+f' index, weight, kx ky kz (crystal coords)\n'
        n_print = min(self.num_k,51)
        for kk in range(n_print):
            message = message+f'   {kk} {self.weights[kk]:2.3f}   {self.kpts[kk,0]: 2.4f}' \
                    f' {self.kpts[kk,1]: 2.4f} {self.kpts[kk,2]: 2.4f}\n'
        print_stdout(message,msg_type='KPOINTS')

    # ------------------------------------------------------------------------------------------------

    def _set_full_grid(self):

        # full grid from -1/2 to 1/2
        self.num_k = self.kpt_grid[0]*self.kpt_grid[1]*self.kpt_grid[2]
        self.kpts = np.zeros((self.num_k,3))
        self.weights = np.ones(self.num_k)
        kx_arr = np.linspace(-1/2,1/2,self.kpt_grid[0])
        ky_arr = np.linspace(-1/2,1/2,self.kpt_grid[1])
        kz_arr = np.linspace(-1/2,1/2,self.kpt_grid[2])
        ind = 0
        for kx in kx_arr:
            for ky in ky_arr:
                for kz in kz_arr:
                    self.kpts[ind,:] = [kx,ky,kz]
                    ind = ind+1

    # ------------------------------------------------------------------------------------------------


