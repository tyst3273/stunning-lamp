import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = r"\usepackage{bm}"


def solve_dynmat(q,g=1,G=1,m=1,M=1,a=1):

    nq = q.size
    freqs = np.zeros((nq,2))
    evecs = np.zeros((nq,2,2),dtype=complex) # nq, nbasis, nmodes

    for ii in range(nq):

        _q = 2*np.pi/a*q[ii]

        D = np.zeros((2,2),dtype=complex)

        D[0,0] = (G+g)/m
        D[0,1] = -(G+g*np.exp(-1j*_q*a))/np.sqrt(m*M)
        D[1,1] = (G+g)/M
        D[1,0] = np.conj(D[0,1])

        _evals, _evecs = eigh(D)
        _freqs = np.sqrt(_evals)

        freqs[ii,:] = _freqs
        evecs[ii,:,:] = _evecs

    return freqs, evecs

N = 81
nq = 101
q = np.linspace(-3/2,3/2,nq)

fig, ax = plt.subplots(1,3,figsize=(8,3),gridspec_kw={'wspace':0.1})

freqs_1, _ = solve_dynmat(q,g=1,G=1,m=1,M=1)
freqs_1 /= 2 #np.sqrt((1+1)*(1+1)/(1*1))

ax[0].plot(q,freqs_1[:,0],marker='o',ms=2,c='r',lw=1,ls='-')
ax[0].plot(q,freqs_1[:,1],marker='o',ms=2,c='b',lw=1,ls='-')

freqs_2, _ = solve_dynmat(q,g=1,G=1.25,m=1,M=1)
freqs_2 /= 2 #np.sqrt((1+1)*(1+1)/(1*1))

ax[1].plot(q,freqs_2[:,0],marker='o',ms=2,c='r',lw=1,ls='-')
ax[1].plot(q,freqs_2[:,1],marker='o',ms=2,c='b',lw=1,ls='-')

freqs_3, _ = solve_dynmat(q,g=1,G=1.5,m=1,M=1)
freqs_3 /= 2 #np.sqrt((1+1)*(1+1)/(1*1))

ax[2].plot(q,freqs_3[:,0],marker='o',ms=2,c='r',lw=1,ls='-')
ax[2].plot(q,freqs_3[:,1],marker='o',ms=2,c='b',lw=1,ls='-')

ax[0].annotate('G=1.00',xycoords='data',xy=(-1.2,1.15),fontsize='large')
ax[1].annotate('G=1.25',xycoords='data',xy=(-1.2,1.15),fontsize='large')
ax[2].annotate('G=1.50',xycoords='data',xy=(-1.2,1.15),fontsize='large')

prim_freqs = 2*np.abs(np.sin(q*np.pi))
prim_freqs /= 2
ax[0].plot(q*2,prim_freqs,marker='o',ms=0,c='k',lw=2,ls=(0,(4,1,2,1)))
ax[1].plot(q*2,prim_freqs,marker='o',ms=0,c='k',lw=2,ls=(0,(4,1,2,1)))
ax[2].plot(q*2,prim_freqs,marker='o',ms=0,c='k',lw=2,ls=(0,(4,1,2,1)))


for ii in range(3):

    ax[ii].plot([-1/2,-1/2],[0,5],lw=1,ls='--',c=(0.2,0.2,0.2))
    ax[ii].plot([1/2,1/2],[0,5],lw=1,ls='--',c=(0.2,0.2,0.2))

    ax[ii].plot([-1,-1],[0,5],lw=1,ls=':',c=(0.2,0.2,0.2))
    ax[ii].plot([1,1],[0,5],lw=1,ls=':',c=(0.2,0.2,0.2))


    for axis in ['top','bottom','left','right']:
        ax[ii].spines[axis].set_linewidth(1.5)

        ax[ii].minorticks_on()
        ax[ii].tick_params(which='both',width=1,labelsize='x-large')
        ax[ii].tick_params(which='major',length=5)
        ax[ii].tick_params(which='minor',length=2)
        #ax[ii].axis([-1/2,1/2,0,2.5])
        ax[ii].set_ylim(0,1.25)
        ax[ii].set_xlim(-1.5,1.5)

        ax[ii].set_xlabel(r'k [$2\pi/a$]',fontsize='x-large')

ax[1].set_yticklabels([])
ax[2].set_yticklabels([])
#ax[2].set_xlabels([])

ax[0].set_ylabel(r'$\omega/\sqrt{\frac{(G+g)(m_1+m_2)}{m_1m_2}}$',fontsize='x-large')



plt.savefig('1d_diatomic_dispersion.png',dpi=300,bbox_inches='tight')


"""
for ii, nn in enumerate(n):

    ax[0].plot(q,freqs[0],marker='o',ms=4,c='r',lw=1,ls='-')
    ax[0].plot(q,freqs[1],marker='o',ms=4,c='b',lw=1,ls='-')

    q = 2*np.pi*nn/N
    cos = np.cos(q*x)
    u = u_q*cos
    _s = np.ones(N)*ii
#    plt.plot(x,_s,marker='o',ms=6,mfc='w',mec='k',mew=1,lw=0)
    plt.plot(x+u,_s,marker='o',ms=8,c='b',mew=0,lw=0)
    _cos = np.cos(q*xf)/5
    plt.plot(xf,_cos+ii,marker='o',ms=0,c='k',mew=0,lw=1,ls='-')
    plt.quiver(x,_s,u,np.zeros(N),angles='xy', scale_units='xy', scale=1,units='xy',
        headwidth=2,headlength=3,headaxislength=3)
    plt.annotate(rf'q=$\frac{{2\pi}}{{a}} \frac{{{nn}}}{{{N}}}$',xycoords='data',
            xy=(7.5,ii-0.2),fontsize='large')
"""

#plt.anotate(

"""
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)

    ax.minorticks_on()
    ax.tick_params(which='both',width=1,labelsize='x-large')
    ax.tick_params(which='major',length=5)
    ax.tick_params(which='minor',length=2)


ax.set_xlabel('i [unitcell label]',fontsize='x-large')
ax.set_yticklabels([])


plt.savefig('1d_harmonic_disp.png',dpi=300,bbox_inches='tight')

"""

#plt.show()

