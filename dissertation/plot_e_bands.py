import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl
from scipy.linalg import eigh

mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = r"\usepackage{bm}"

a = 1

nk = 251
k = np.linspace(-1/2,1/2,nk,endpoint=False)

nx = 100
x = np.linspace(0,a,nx,endpoint=False)

n = 1
V0 = 5
V = V0*np.cos(2*np.pi*n*x)/2 

VG = np.fft.fftfreq(nx,x[1]-x[0])
V_fft = np.fft.fft(V,norm='forward').round(3)

Gcut = VG.max()/2
G = VG[np.flatnonzero(np.abs(VG) <= Gcut)] 
nG = G.size

KE = np.zeros((nG,nG))
PE = np.zeros((nG,nG),dtype=complex)

nb = 2
ek = np.zeros((nk,nb))
wk = np.zeros((nk,nb,nG),dtype=complex)

ek0 = np.zeros((nk,nb))
wk0 = np.zeros((nk,nb,nG),dtype=complex)

for kk in range(nk):

    KE[...] = 0.0
    PE[...] = 0.0

    for ii in range(nG):
        for jj in range(ii,nG):
            
            if ii == jj:
                KE[ii,jj] = 4*np.pi**2*0.5*(k[kk]+G[ii])**2
            
            Gii = G[ii]; Gjj = G[jj]

            GG = Gii-Gjj
            ind = np.flatnonzero(VG == GG)[0]
            PE[ii,jj] = V_fft[ind]

    evals, evecs = eigh(KE,subset_by_index=[0,nb-1])
    ek0[kk,:] = evals
    wk0[kk,...] = evecs.T

    evals, evecs = eigh(KE+PE,subset_by_index=[0,nb-1],lower=False)
    ek[kk,:] = evals
    wk[kk,...] = evecs.T


fig, ax = plt.subplots(1,2,figsize=(9,4),gridspec_kw={'wspace':0.25})

ax[0].plot(x,np.zeros(nx),c='r',lw=1,ls=(0,(4,2,2,2)),label=r'$V=0$')
ax[0].plot(x,V,c='b',lw=1,ms=0,label=r'$V=V_0 \cos(2\pi x)$')

ax[1].plot(k,ek0/V0,c='r',lw=1,ls=(0,(4,2,2,2)),ms=0)
ax[1].plot(k,ek/V0,c='b',lw=1,ms=0)

for ii in range(2):

    for axis in ['top','bottom','left','right']:
        ax[ii].spines[axis].set_linewidth(1.1)

    ax[ii].minorticks_on()
    ax[ii].tick_params(which='both', width=1, labelsize='x-large')
    ax[ii].tick_params(which='major', length=5)
    ax[ii].tick_params(which='minor', length=3, color='k')
    #ax[ii].set_rasterization_zorder(10000)

ax[0].set_xlabel(r'$x/a$',fontsize='large')
ax[0].set_ylabel(r'$V/V_0$',labelpad=3.0,fontsize='large')

ax[1].set_xlabel(r'$k~[2\pi/a]$',fontsize='large')
ax[1].set_ylabel(r'$\epsilon/V_0$',labelpad=10.0,fontsize='large')

ax[0].annotate('(a)',fontsize='large',xy=(0.01,0.025),xycoords='axes fraction')
ax[1].annotate('(b)',fontsize='large',xy=(0.01,0.025),xycoords='axes fraction')

ax[0].legend(frameon=False,loc='center',bbox_to_anchor=(0.5,0.75))

#ax.autoscale(tight=True)
#ax.axis('tight')

#ax.axis([0.1,1000,0.1,1000])

plt.savefig('e_bands_model.pdf',dpi=100,bbox_inches='tight')
plt.show()
