import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl
from solve_1d_e_bands import solve

mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = r"\usepackage{bm}"

ek0, wk0, ek, wk, x, k, G, V, V0 = solve()
nk = k.size
nx = x.size
nG = G.size
nb = ek0.shape[1]

fig, ax = plt.subplots(1,2,figsize=(9,4),gridspec_kw={'wspace':0.25})

ax[0].plot(x,np.zeros(nx),c='r',lw=1,ls=(0,(4,2,2,2)),label=r'$V=0$')
ax[0].plot(x,V/V0,c='b',lw=1,ms=0,label=r'$V=V_0 \cos(2\pi x)$')

ax[1].plot(k,ek0/V0,c='r',lw=1,ls=(0,(4,2,2,2)),ms=0)
ax[1].plot(k,ek/V0,c='b',lw=1,ms=0)

for ii in range(2):

    for axis in ['top','bottom','left','right']:
        ax[ii].spines[axis].set_linewidth(1.1)

    ax[ii].minorticks_on()
    ax[ii].tick_params(which='both', width=1, labelsize='x-large')
    ax[ii].tick_params(which='major', length=5)
    ax[ii].tick_params(which='minor', length=3, color='k')

ax[0].set_xlabel(r'$x/a$',fontsize='large')
ax[0].set_ylabel(r'$V/V_0$',labelpad=3.0,fontsize='large')

ax[1].set_xlabel(r'$k~[2\pi/a]$',fontsize='large')
ax[1].set_ylabel(r'$\epsilon/V_0$',labelpad=10.0,fontsize='large')

ax[0].annotate('(a)',fontsize='large',xy=(0.01,0.025),xycoords='axes fraction')
ax[1].annotate('(b)',fontsize='large',xy=(0.01,0.025),xycoords='axes fraction')

ax[0].legend(frameon=False,loc='center',bbox_to_anchor=(0.5,0.75))

ax[0].axis([0,1,-0.5,0.5])
ax[1].axis([-0.5,0.5,-0.1,4])

plt.savefig('e_bands_model.pdf',dpi=100,bbox_inches='tight')
#plt.show()
