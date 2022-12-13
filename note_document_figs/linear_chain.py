

import numpy as np
import matplotlib.pyplot as plt



def get_omega(q,K,m):
    return 2*np.sqrt(K/m)*np.abs(np.sin(q/2))

def get_omega_diatomic(q,G,g,M,m):
    _p = (G+g)*(m+M)/(2*m*M)+np.sqrt((G+g)**2*(m+M)**2-16*G*g*M*m*np.sin(q/2)**2)/(2*m*M)
    _m = (G+g)*(m+M)/(2*m*M)-np.sqrt((G+g)**2*(m+M)**2-16*G*g*M*m*np.sin(q/2)**2)/(2*m*M)
    return np.sqrt(_p), np.sqrt(_m)


nq = 51
q = np.linspace(-1,1,nq)*np.pi

m = 1
M = 1.5*m

K = 2

g = 1
G = 1.5*g

e1d = np.zeros(nq)
e2d = np.zeros((nq,2))

for _iq, _q in enumerate(q):
    e1d[_iq] = get_omega(_q,K,m)
    e2d[_iq,:] = get_omega_diatomic(_q,G,g,M,m)

q = np.linspace(-1,1,nq)
plt.plot(q,e1d/np.sqrt(K/m),c='k',marker='o',ms=4,mfc='k',mec='k',lw=1,ls='-')
plt.xlabel(r'q [2$\pi$/a]')
plt.ylabel(r'$\omega$/$\sqrt{K/m}$')
#plt.show()
plt.savefig('1d_harmonic_chain.png',bbox_inches='tight',dpi=150)
plt.close()

f = np.sqrt((G+g)*(m+M)/(m*M))
plt.plot(q,e2d[:,0]/f,c='b',marker='o',ms=4,mfc='b',mec='b',lw=1,ls='-')
plt.plot(q,e2d[:,1]/f,c='r',marker='o',ms=4,mfc='r',mec='r',lw=1,ls='-')
plt.xlabel(r'q [2$\pi$/a]')
plt.ylabel(r'$\omega$/$\sqrt{(G+g)(m+M)/mM}$')
#plt.show()
plt.savefig('1d_diatomic_chain.png',bbox_inches='tight',dpi=150)
plt.close()




