import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = r"\usepackage{bm}"

a = 1

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

m = 1
M = 2
g = 1
G = 2

u_q = 0.25

N = 8
n = [0,1,4]
q = np.array(n)/N

fig, ax = plt.subplots(figsize=(8,3))
freqs, evecs = solve_dynmat(q,g=g,G=G,m=m,M=M)

for ii in range(len(n)):
    for jj in range(2):
        
        e = evecs[ii,:,jj] 
        ind = np.argmax(np.abs(e))
        phase = e[ind]/np.abs(e[ind])
        e /= phase
        e = np.real(e)/e.max()
        evecs[ii,:,jj] = e

evecs = np.real(evecs)

evecs[:,0,:] /= np.sqrt(m)
evecs[:,1,:] /= np.sqrt(M)

u_q = evecs*0.25

x1 = np.arange(N)
x2 = np.arange(N)+1/2
xf = np.linspace(0,N,101)

q = 2*np.pi*q/a

for ii in range(N):
    plt.plot([ii,ii],[0-0.5,len(n)-0.5],lw=1,ls='--',c=(0.2,0.2,0.2))
    plt.plot([ii+1/2,ii+1/2],[0-0.5,len(n)-0.5],lw=1,ls=':',c=(0.2,0.2,0.2))

for ii in range(len(n)):

    cos = np.exp(1j*q[ii]*x1)
    
    u11 = np.real(u_q[ii,0,0]*cos)
    u21 = np.real(u_q[ii,1,0]*cos)

    u12 = np.real(u_q[ii,0,1]*cos)
    u22 = np.real(u_q[ii,1,1]*cos)

    print(u12)

    _s = np.ones(N)*ii
    plt.plot(x1,_s,marker='o',ms=6,mfc='w',mec='k',mew=1,lw=0)

    plt.plot(x1+u11,_s,marker='o',ms=6,c='m',mew=0,lw=0)
    plt.plot(x2+u21,_s,marker='o',ms=8,c='r',mew=0,lw=0)

    plt.plot(x1+u12,_s+0.25,marker='o',ms=6,c='m',mew=0,lw=0)
    plt.plot(x2+u22,_s+0.25,marker='o',ms=8,c='r',mew=0,lw=0)


    _cos = np.cos(q[ii]*xf)/5
    plt.plot(xf,_cos+ii,marker='o',ms=0,c='k',mew=0,lw=1,ls='-')

    plt.quiver(x1,_s,u11,np.zeros(N),angles='xy', scale_units='xy', scale=1,units='xy',
        headwidth=2,headlength=3,headaxislength=3)
    plt.quiver(x2,_s,u21,np.zeros(N),angles='xy', scale_units='xy', scale=1,units='xy',
        headwidth=2,headlength=3,headaxislength=3)
    plt.quiver(x1,_s+0.25,u12,np.zeros(N),angles='xy', scale_units='xy', scale=1,units='xy',
        headwidth=2,headlength=3,headaxislength=3)
    plt.quiver(x2,_s+0.25,u22,np.zeros(N),angles='xy', scale_units='xy', scale=1,units='xy',
        headwidth=2,headlength=3,headaxislength=3)

    plt.annotate(rf'q=$\frac{{2\pi}}{{a}} \frac{{{n[ii]}}}{{{N}}}$',xycoords='data',
            xy=(7.6,ii-0.25),fontsize='large')

#plt.anotate(

for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)

    ax.minorticks_on()
    ax.tick_params(which='both',width=1,labelsize='x-large')
    ax.tick_params(which='major',length=5)
    ax.tick_params(which='minor',length=2)


ax.set_xlabel('i [unitcell label]',fontsize='x-large')
ax.set_yticklabels([])


plt.savefig('1d_diatomic_displacement.png',dpi=300,bbox_inches='tight')


plt.show()

