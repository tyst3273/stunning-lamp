
import numpy as np
import matplotlib.pyplot as plt

A = 0.4582
B = 0.44
C = 7.8

def rs(n):

    return (3/(4*np.pi*n))**(1/3)

def ex(r):

    return  - A / r


def ec(r):

    return - B / ( r + C )

def deriv(r):

    return - 4 * np.pi / 9 * ( 0.4582 * r**2 + 0.44 * r**4 / (r+7.8)**2 )

fig, ax = plt.subplots(1,2)



x = np.linspace(-1,2,1000)

ai = 10
aj = 10

xi = 0
xj = 1

ni = (ai/np.pi)**(3/2) * np.exp(-ai * (x-xi)**2 )
nj = (aj/np.pi)**(3/2) * np.exp(-aj * (x-xj)**2 )
n = ni+nj

aij = ai+aj
aiaj = ai*aj
rij = (ai * xi + aj*xj) / aij

nij = (ai*aj/np.pi**2)**(3/2) * np.exp(-aiaj * (xi-xj)**2 / aij) * np.exp( -aij * ( x - rij )**2 )

ax[0].axvline(rij,lw=0.5,ls=(0,(2,1)),c='k')

ax[0].plot(x,n / 3,c='r',lw=1,label=r'$\rho_I+\rho_J$')
ax[0].plot(x,ni / 3,c='r',lw=1,ls=(0,(2,1)),label=r'$\rho_{I,J}$')
ax[0].plot(x,nj / 3,c='r',lw=1,ls=(0,(2,1)))

ax[0].plot(x,nij / 3**2,c='m',lw=1,label=r'$\rho_I\rho_J$')
# ax[0].plot(n_m**2+n_p**2,c='m',lw=1,ls=(0,(2,1)))

r = rs(n)
K = deriv(r)

ax[0].plot(x,K,c='k',lw=1,label=r'$K(\rho_I+\rho_J)$')
ax[0].plot(x,K*nij,c='b',lw=1,label=r'$K(\rho_I+\rho_J) \rho_I \rho_J$')


ax[0].set_title(r'$\alpha_I=$'+f'{ai:.1f}, '+r'$\alpha_J=$'+f'{aj:.1f}')



ai = 10
aj = 5

xi = 0
xj = 1

ni = (ai/np.pi)**(3/2) * np.exp(-ai * (x-xi)**2 )
nj = (aj/np.pi)**(3/2) * np.exp(-aj * (x-xj)**2 )
n = ni+nj

aij = ai+aj
aiaj = ai*aj
rij = (ai * xi + aj*xj) / aij

nij = (ai*aj/np.pi**2)**(3/2) * np.exp(-aiaj * (xi-xj)**2 / aij) * np.exp( -aij * ( x - rij )**2 )

ax[1].axvline(rij,lw=0.5,ls=(0,(2,1)),c='k')

ax[1].plot(x,n / 6,c='r',lw=1)
ax[1].plot(x,ni / 3,c='r',lw=1,ls=(0,(2,1)))
ax[1].plot(x,nj / 3,c='r',lw=1,ls=(0,(2,1)))

ax[1].plot(x,nij / 3**2,c='m',lw=1)
# ax[1].plot(n_m**2+n_p**2,c='m',lw=1,ls=(0,(2,1)))

r = rs(n)
K = deriv(r)

ax[1].plot(x,K,c='k',lw=1)
ax[1].plot(x,K*nij,c='b',lw=1)

ax[1].set_title(r'$\alpha_I=$'+f'{ai:.1f}, '+r'$\alpha_J=$'+f'{aj:.1f}')


fig.legend(ncol=5,loc='upper center',bbox_to_anchor=[0.5,1.025])
ax[0].axis([-1,2,-0.5,0.5])
ax[1].axis([-1,2,-0.5,0.5])

ax[0].set_xlabel('|r|')
ax[1].set_xlabel('|r|')

plt.savefig('E2C_advanced.png',format='png',dpi=300)

plt.show()