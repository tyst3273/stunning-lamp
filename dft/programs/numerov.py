import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------------------------------------

class numerov_integrator:

    # -----------------------------------------------------------------------------------------------------

    def __init__(self):
        
        self.const = 1 # 2m/hbar**2

    # -----------------------------------------------------------------------------------------------------

    def init_grid(self,x_max,x_min=0,dx=0.01):

        self.x_min = x_min
        self.x_max = x_max
        self.dx = dx
        self.x_grid = np.arange(x_min,x_max+dx,dx)
        self.num_grid = self.x_grid.size

    # -----------------------------------------------------------------------------------------------------

    def numerov(self,pot,x0=0,y0=2,y1=1,xf=8,e_val=1):

        """
        integrates to the right
        """

        f = 1+self.const*(e_val-pot)*self.dx**2/12

        self.sol = np.zeros(self.num_grid)

        n = np.argwhere(self.x_grid >= x0).flatten().min()
        nf = np.argwhere(self.x_grid <= xf).flatten().max()

        n_1 = n+1
        n_2 = n_1+1

        self.start = n 
        self.end = nf

        self.sol[n] = y0
        self.sol[n_1] = y1

        num_to_do = nf-n-2
        print(num_to_do)

        for ii in range(num_to_do):

            self.sol[n_2] = ((12-10*f[n_1])*self.sol[n_1]-f[n]*self.sol[n])/f[n_2]
            n = n_1
            n_1 = n_2
            n_2 = n_2+1


    # -----------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------
# run the thing

if __name__ == '__main__':

    numerov = numerov_integrator()
    numerov.init_grid(x_max=10,x_min=-10,dx=0.001)

    k = 1
    pot = k/2*numerov.x_grid**2

    numerov.numerov(pot,x0=0,y0=20,y1=15,e_val=10,xf=8)

    plt.plot(numerov.x_grid,pot,color='b')
    plt.plot(numerov.x_grid,numerov.sol,color='r')
    plt.show()





















