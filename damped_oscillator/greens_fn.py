
import numpy as np
import matplotlib.pyplot as plt


class c_greens_fn:

    # ----------------------------------------------------------------------------------------------

    def __init__(self,t=None,m=1,w0=1,gamma=0):

        self.m = m
        self.gamma = gamma
        self.w0 = w0

        if t is None:
            self.t = np.linspace(0,10,1001)
        else:
            self.t = np.array(t)

        self._get_green_fn()

    # ----------------------------------------------------------------------------------------------

    def _get_green_fn(self):
        
        gamma = self.gamma
        w0 = self.w0
        t = self.t

        print(gamma,w0)

        damp = np.exp(-gamma*t)

        if gamma < w0:
            c = np.sqrt(w0**2-gamma**2)
            g = np.sin(c*t)/c
        elif gamma > w0:
            c = np.sqrt(gamma**2-w0**2)
            g = np.sinh(c*t)/c
        else:
            g = t

        mask = (t >= 0).astype(float)
        self.greens_fn = g*mask*damp

    # ----------------------------------------------------------------------------------------------

    def plot_greens_fn(self):

        plt.plot(t,self.greens_fn)
        plt.xlabel('t')
        plt.ylabel('G(t)')
        plt.show()
    
    # ----------------------------------------------------------------------------------------------


if __name__ == '__main__':

    m = 1
    gamma = 0.01
    w0 = 5

    t_max = 7*np.pi/w0
    num_t = 1001
    t = np.linspace(0,t_max,num_t)

    g = c_greens_fn(t=t,m=m,gamma=gamma,w0=w0)
    g = g.greens_fn

    w = np.linspace(0,3*w0,num_t)

    a = 2*gamma*w/((w**2-w0**2)**2+(2*gamma*w)**2)

    hi = max([g.max(),a.max()])
    lo = min([g.min(),a.min()])

    a *= hi/a.max()
    g *= hi/g.max()

    fig, ax = plt.subplots()
    ax2 = ax.twiny()

    ax.plot(w,a,c='b',label=r'$G(\omega)$')
    ax.plot([w0,w0],[lo,hi],ls='--',c='k')

    ax2.plot(t,g,c='r',label='$G(t)$')

    ax.set_ylabel('Im(G) [arb units]')
    ax2.set_xlabel(r'$t$')
    ax.set_xlabel(r'$\omega$')

    fig.legend()

    ax.set_title(r'$\omega_0$='+f'{w0:3.2f},  '+r'$\gamma$='+f'{gamma:3.2f}')

    plt.savefig(f'w0_{w0}_gamma_{gamma}.pdf',dpi=150,bbox_inches='tight')
    plt.show()
    








