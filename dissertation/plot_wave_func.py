import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl
from scipy.linalg import eigh

mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = r"\usepackage{bm}"


nk = 251
k = np.linspace(-1/2,1/2,nk,endpoint=False)

nx = 51
x = np.linspace(0,1,nx,endpoint=False)

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
                KE[ii,jj] = 4*np.pi**2*0.5*(k[kk]+G[ii])
            
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


nxp = 501
xp = np.linspace(0,2,nxp,endpoint=False)
#Vp = V0*np.cos(2*np.pi*n*xp)/2

u = np.zeros((nk,nb,nxp),dtype=complex)
u0 = np.zeros((nk,nb,nxp),dtype=complex)

for kk in range(nk):
    for ii in range(nb):
        for jj in range(nG):        
            
            u[kk,ii,:] += wk[kk,ii,jj]*np.exp(2j*np.pi*G[jj]*xp) 
            u0[kk,ii,:] += wk0[kk,ii,jj]*np.exp(2j*np.pi*G[jj]*xp)

rho = np.sum(np.abs(u)**2,axis=0)/nk
rho0 = np.sum(np.abs(u0)**2,axis=0)/nk

#ax.plot(xp,np.zeros(xp.size),c='k',lw=1,ls=(0,(4,2,2,2)),label=r'$V=0$')
#ax.plot(xp,Vp,c='m',lw=1,ms=0,label=r'$V=V_0 \cos(2\pi x)$')

fig, ax = plt.subplots(figsize=(9,4),gridspec_kw={'wspace':0.25})

#ax.plot(xp,rho[0,:],c='r',lw=1,ms=0,ls=(0,(4,2,2,2)))
#ax.plot(xp,rho0[0,:],c='b',lw=1,ms=0)
#ax.fill_between(xp,0,rho[0,:],color='r',alpha=0.1)
#ax.fill_between(xp,0,rho0[0,:],color='b',alpha=0.1)

ind = np.flatnonzero(k == -0.5)[0]
uu = u[ind,0,:]
phase = uu[0]/np.abs(uu[0])
uu /= phase
w = uu*np.exp(2j*np.pi*k[0]*xp) #/nk
plt.plot(xp,np.real(uu),c='b',lw=1,ms=0)
plt.plot(xp,np.exp(2j*np.pi*k[0]*xp),c='k',lw=1)
plt.plot(xp,np.real(w),c='m',lw=1,ms=0)


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

#ax.axis([0.1,1000,0.1,1000])

plt.savefig('e_bands_wavefunctions.pdf',dpi=100,bbox_inches='tight')
plt.show()
