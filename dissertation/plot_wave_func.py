import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl
from scipy.linalg import eigh
from solve_1d_e_bands import solve

mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = r"\usepackage{bm}"

ek0, wk0, ek, wk, x, k, G, V, V0 = solve(nb=25)
nk = k.size
nx = x.size
nG = G.size
nb = ek0.shape[1]

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

fig, ax = plt.subplots(1,2,figsize=(9,4),gridspec_kw={'wspace':0.1})

# -------------------
# no PE

kpt = -0.5
k_ind = np.flatnonzero(k == kpt)[0]
exp = np.exp(2j*np.pi*k[k_ind]*xp)

ax[0].plot(xp,np.zeros(nxp),c='m',lw=1,ms=0) #label=r'$V=V_0 \cos(2\pi x)$')
ax[0].plot(xp,np.real(exp),c='m',lw=1,ls=(0,(4,1,2,1)))
    
for ii in range(nc):
    ax[0].plot([ii,ii],[-10,10],ls=(0,(1,1)),c=(0.5,0.5,0.5),ms=0,lw=1)

s = 4
ax[0].plot(xp,rho0[0,:]+s,c='r',lw=1,ms=0,ls=(0,(4,2,2,2)))
ax[0].fill_between(xp,s,rho0[0,:]+s,color='r',alpha=0.1)
ax[0].plot([0,nc],[s,s],lw=1,ls=(0,(1,1)),ms=0,c=(0.5,0.5,0.5))

uu = u0[k_ind,0,:]
ind = np.argmax(np.abs(uu))
phase = uu[ind]/np.abs(uu[ind])
uu /= phase
w = uu*exp #/nk
ind = np.argmax(np.abs(w))
phase = w[ind]/np.abs(w[ind])
w /= phase
ax[0].plot(xp,np.real(uu)+s,c='r',lw=1,ms=0)
ax[0].plot(xp,np.real(w)+s,c='r',lw=1,ms=0,ls=(0,(2,1)))

s = 7
ax[0].plot(xp,rho0[1,:]+s,c='b',lw=1,ms=0,ls=(0,(4,1,2,1)))
ax[0].fill_between(xp,s,rho0[1,:]+s,color='b',alpha=0.1)
ax[0].plot([0,nc],[s,s],lw=1,ls=(0,(1,1)),ms=0,c=(0.5,0.5,0.5))

uu = u0[k_ind,1,:]
ind = np.argmax(np.abs(uu))
phase = uu[ind]/np.abs(uu[ind])
uu /= phase
w = uu*exp #/nk
ind = np.argmax(np.abs(w))
phase = w[ind]/np.abs(w[ind])
w /= phase
ax[0].plot(xp,np.real(uu)+s,c='b',lw=1,ms=0)
ax[0].plot(xp,np.real(w)+s,c='b',lw=1,ms=0,ls=(0,(2,1)))

s = 10
ax[0].plot(xp,rho0[11,:]+s,c='g',lw=1,ms=0,ls=(0,(4,1,2,1)))
ax[0].fill_between(xp,s,rho0[11,:]+s,color='g',alpha=0.1)
ax[0].plot([0,nc],[s,s],lw=1,ls=(0,(1,1)),ms=0,c=(0.5,0.5,0.5))

uu = u0[k_ind,11,:]
ind = np.argmax(np.abs(uu))
phase = uu[ind]/np.abs(uu[ind])
uu /= phase
w = uu*exp #/nk
ind = np.argmax(np.abs(w))
phase = w[ind]/np.abs(w[ind])
w /= phase
ax[0].plot(xp,np.real(uu)+s,c='g',lw=1,ms=0)
ax[0].plot(xp,np.real(w)+s,c='g',lw=1,ms=0,ls=(0,(2,1)))

# -------------------
# PE

kpt = -0.5
k_ind = np.flatnonzero(k == kpt)[0]
exp = np.exp(2j*np.pi*k[k_ind]*xp)

ax[1].plot(xp,Vp/V0,c='m',lw=1,ms=0) #,label=r'$V=V_0 \cos(2\pi x)$')
ax[1].plot(xp,np.real(exp),c='m',lw=1,ls=(0,(4,1,2,1)))

for ii in range(nc):
    ax[1].plot([ii,ii],[-10,10],ls=(0,(1,1)),c=(0.5,0.5,0.5),ms=0,lw=1)

s = 4
ax[1].plot(xp,rho[0,:]+s,c='r',lw=1,ms=0,ls=(0,(4,1,2,1)))
ax[1].fill_between(xp,s,rho[0,:]+s,color='r',alpha=0.1)
ax[1].plot([0,nc],[s,s],lw=1,ls=(0,(1,1)),ms=0,c=(0.5,0.5,0.5))

uu = u[k_ind,0,:]
ind = np.argmax(np.abs(uu))
phase = uu[ind]/np.abs(uu[ind])
uu /= phase
w = uu*exp #/nk
ind = np.argmax(np.abs(w))
phase = w[ind]/np.abs(w[ind])
w /= phase
ax[1].plot(xp,np.real(uu)+s,c='r',lw=1,ms=0)
ax[1].plot(xp,np.real(w)+s,c='r',lw=1,ms=0,ls=(0,(2,1)))

s = 7
ax[1].plot(xp,rho[1,:]+s,c='b',lw=1,ms=0,ls=(0,(4,1,2,1)))
ax[1].fill_between(xp,s,rho[1,:]+s,color='b',alpha=0.1)
ax[1].plot([0,nc],[s,s],lw=1,ls=(0,(1,1)),ms=0,c=(0.5,0.5,0.5))

uu = u[k_ind,1,:]
ind = np.argmax(np.abs(uu))
phase = uu[ind]/np.abs(uu[ind])
uu /= phase
w = uu*exp #/nk
ind = np.argmax(np.abs(w))
phase = w[ind]/np.abs(w[ind])
w /= phase
ax[1].plot(xp,np.real(uu)+s,c='b',lw=1,ms=0)
ax[1].plot(xp,np.real(w)+s,c='b',lw=1,ms=0,ls=(0,(2,1)))

s = 10
ax[1].plot(xp,rho[11,:]+s,c='g',lw=1,ms=0,ls=(0,(4,1,2,1)))
ax[1].fill_between(xp,s,rho[11,:]+s,color='g',alpha=0.1)
ax[1].plot([0,nc],[s,s],lw=1,ls=(0,(1,1)),ms=0,c=(0.5,0.5,0.5))

uu = u[k_ind,11,:]
ind = np.argmax(np.abs(uu))
phase = uu[ind]/np.abs(uu[ind])
uu /= phase
w = uu*exp #/nk
ind = np.argmax(np.abs(w))
phase = w[ind]/np.abs(w[ind])
w /= phase
ax[1].plot(xp,np.real(uu)+s,c='g',lw=1,ms=0)
ax[1].plot(xp,np.real(w)+s,c='g',lw=1,ms=0,ls=(0,(2,1)))


# ------------
# format plot

l1, = ax[0].plot([-2,-1],[0,0],lw=1,ms=0,ls=(0,(4,1,2,1)),c='g')
l2, = ax[0].plot([-2,-1],[0,0],lw=1,ms=0,ls=(0,(2,1)),c='g')
l3, = ax[0].plot([-2,-1],[0,0],lw=1,ms=0,ls='-',c='g')
labels = [r'$u_{kn}(x)$',r'$\psi_{kn}(x)$',r'$\rho_{n}(x)$']

leg1 = ax[0].legend([l1,l2,l3],labels,frameon=False,loc='upper center',bbox_to_anchor=(0.5,1),ncols=3,
    labelspacing=0.1,handletextpad=0.1,columnspacing=0.5,handlelength=2)

l1, = ax[1].plot([-2,-1],[0,0],lw=1,ms=0,ls=(0,(4,1,2,1)),c='m')
l2, = ax[1].plot([-2,-1],[0,0],lw=1,ms=0,ls='-',c='m')
labels = [r'$\exp(ikx)$',r'$V(x)/V_0$']

ax[0].legend([l1,l2],labels,frameon=False,loc='upper center',bbox_to_anchor=(0.5,0.275),ncols=3,
    labelspacing=0.1,handletextpad=0.1,columnspacing=0.5,handlelength=2)
ax[0].add_artist(leg1)

fig.suptitle(r'$k=1/2$',fontsize='x-large',y=0.98)

ax[1].annotate(r'n=1',fontsize='large',xy=(0.45,0.275),xycoords='axes fraction')
ax[1].annotate(r'n=2',fontsize='large',xy=(0.45,0.5),xycoords='axes fraction')
ax[1].annotate(r'n=12',fontsize='large',xy=(0.45,0.9),xycoords='axes fraction')


for ii in range(2):

    for axis in ['top','bottom','left','right']:
        ax[ii].spines[axis].set_linewidth(1.1)

    ax[ii].minorticks_on()
    ax[ii].tick_params(which='both', width=1, labelsize='x-large')
    ax[ii].tick_params(which='major', length=5)
    ax[ii].tick_params(which='minor', length=3, color='k')

ax[0].set_xlabel(r'$x/a$',fontsize='large')
ax[1].set_xlabel(r'$x/a$',fontsize='large')

ax[0].set_yticklabels([])
ax[1].set_yticklabels([])

#ax.set_ylabel(r'$V/V_0$',labelpad=3.0,fontsize='large')

#ax[0].legend(frameon=False,loc='center',bbox_to_anchor=(0.5,0.75))

#ax.autoscale(tight=True)
#ax.axis('tight')

ax[0].axis([0,nc,-1.5,13])
ax[1].axis([0,nc,-1.5,13])

plt.savefig(f'e_bands_wavefunctions.pdf',dpi=100,bbox_inches='tight')
#plt.show()
