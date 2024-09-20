import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl
from scipy.linalg import eigh
from solve_1d_e_bands import solve

mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = r"\usepackage{bm}"

ek0, wk0, ek, wk, x, k, G, V, V0 = solve()
nk = k.size
nx = x.size
nG = G.size
nb = ek0.shape[1]

kpt = -0.5
k_ind = np.flatnonzero(k == kpt)[0]

nxp = 501
nc = 4
xp = np.linspace(0,nc,nxp,endpoint=False)
Vp = V0*np.cos(2*np.pi*xp)/2

u = np.zeros((nk,nb,nxp),dtype=complex)
u0 = np.zeros((nk,nb,nxp),dtype=complex)

for kk in range(nk):
    for ii in range(nb):
        for jj in range(nG):        
            
            u[kk,ii,:] += wk[kk,ii,jj]*np.exp(2j*np.pi*G[jj]*xp) 
            u0[kk,ii,:] += wk0[kk,ii,jj]*np.exp(2j*np.pi*G[jj]*xp)

rho = np.sum(np.abs(u)**2,axis=0)/nk
rho0 = np.sum(np.abs(u0)**2,axis=0)/nk

fig, ax = plt.subplots(figsize=(9,4),gridspec_kw={'wspace':0.25})

exp = np.exp(2j*np.pi*k[k_ind]*xp)

#ax.plot(xp,np.zeros(xp.size),c='k',lw=1,ls=(0,(4,2,2,2)),label=r'$V=0$')
ax.plot(xp,Vp/V0,c='k',lw=1,ms=0,label=r'$V=V_0 \cos(2\pi x)$')
plt.plot(xp,np.real(exp),c='m',lw=1,ls=(0,(4,1,2,1)))

for ii in range(nc):
    ax.plot([ii,ii],[-10,10],ls=(0,(2,1)),c=(0.25,0.25,0.25),ms=0,lw=1)

s = 2
ax.plot(xp,rho[0,:]+s,c='r',lw=1,ms=0,ls=(0,(4,2,2,2)))
ax.fill_between(xp,s,rho[0,:]+s,color='r',alpha=0.1)
ax.plot([0,nc],[s,s],lw=1,ls=(0,(2,1)),ms=0,c=(0.25,0.25,0.25))

uu = u[k_ind,0,:]
ind = np.argmax(np.abs(uu))
phase = uu[ind]/np.abs(uu[ind])
uu /= phase
w = uu*exp #/nk
plt.plot(xp,np.real(uu)+s,c='r',lw=1,ms=0)
#plt.plot(xp,np.real(w)+s,c='m',lw=1,ms=0)

s = 6
ax.plot(xp,rho[1,:]+s,c='b',lw=1,ms=0,ls='-')
ax.fill_between(xp,s,rho[1,:]+s,color='b',alpha=0.1)
ax.plot([0,nc],[s,s],lw=1,ls=(0,(2,1)),ms=0,c=(0.25,0.25,0.25))

uu = u[k_ind,1,:]
ind = np.argmax(np.abs(uu))
phase = uu[ind]/np.abs(uu[ind])
uu /= phase
w = uu*exp #/nk
plt.plot(xp,np.real(uu)+s,c='b',lw=1,ms=0)
#plt.plot(xp,np.real(w)+s,c='m',lw=1,ms=0)


for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.1)

ax.minorticks_on()
ax.tick_params(which='both', width=1, labelsize='x-large')
ax.tick_params(which='major', length=5)
ax.tick_params(which='minor', length=3, color='k')
#ax.set_rasterization_zorder(10000)

ax.set_xlabel(r'$x/a$',fontsize='large')
ax.set_ylabel(r'$V/V_0$',labelpad=3.0,fontsize='large')

#ax[0].legend(frameon=False,loc='center',bbox_to_anchor=(0.5,0.75))

#ax.autoscale(tight=True)
#ax.axis('tight')

ax.axis([0,nc,-1,8])

plt.savefig(f'e_bands_wavefunctions_{kpt:3.2f}.pdf',dpi=100,bbox_inches='tight')
#plt.show()
