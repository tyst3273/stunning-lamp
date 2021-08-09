import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------------------------------------

class radial_equation:

    """
    set up radial equation and integrate numerical to get the solution.
    """

    # ----------------------------------------------------------------------------------------------------

    def __init__(self,ang_mom_l,input_energy,radial_grid_step,radial_grid_max,radial_grid_min=0):

        """
        get and set the input options
        """
        
        # set the variables
        self.ang_mom_l = ang_mom_l
        self.input_energy = input_energy
        self.radial_grid_step = radial_grid_step
        self.radial_grid_max = radial_grid_max
        self.radial_grid_min = radial_grid_min # defaults to the origin, r=0.

        # set up the radial grid
        self._set_radial_grid()

    # ----------------------------------------------------------------------------------------------------

    def _set_radial_grid(self):

        """
        set up the grid
        """

        # the radial grid
        self.radial_grid = np.arange(self.radial_grid_min,self.radial_grid_max+self.radial_grid_step,
                              self.radial_grid_step)
        self.num_grid_points = self.radial_grid.shape[0]

        print(self.num_grid_points)

    # ----------------------------------------------------------------------------------------------------

    def numerov(self):

        """
        integrate the radial equation using numerovs method
        """

        pass

    # ----------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------

if __name__ == "__main__":

    ang_mom_l = 0
    input_energy = 1
    radial_grid_step = 0.01
    radial_grid_max = 10

    rad_eq = radial_equation(ang_mom_l=ang_mom_l,input_energy=input_energy,radial_grid_step=radial_grid_step,
                            radial_grid_max=radial_grid_max)





