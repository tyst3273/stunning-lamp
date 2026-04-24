
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

n = np.linspace(0,100,100000)
dn = n[1]-n[0]
r = rs(n)

exc = ex(r)+ec(r)
plt.plot(n,exc,c='m',lw=1,label=r'$\epsilon_{xc}$')

dexc = np.gradient(exc,dn)
plt.plot(n,dexc,c='r',lw=1,label=r'$\partial \epsilon_{xc} / \partial n$')

d2exc = np.gradient(dexc,dn)
plt.plot(n,d2exc,c='b',lw=1,label=r'$\partial^2 \epsilon_{xc} / \partial n^2$')

f = 2*exc + n * d2exc 

plt.plot(n,f,c='k',lw=1,label=r'$\partial^2 \epsilon_{xc} / \partial n^2$')

plt.axhline(0,lw=0.5,ls=(0,(2,1)),c='k')
plt.axvline(0,lw=0.5,ls=(0,(2,1)),c='k')

# I = np.zeros_like(K)

# for x0 in [-2,-1,0,1,2]:
#     n0 = (a/np.pi)**(3/2) * np.exp(-a * (x-x0)**2 )
#     r0 = rs(n0)
#     K0 = ex(r0)+ec(r0)
#     plt.plot(x,K0,c='k',lw=1,ls=(0,(2,1)))
#     plt.plot(x,K0*n0**2,c='b',lw=1,ls=(0,(2,1)))
#     I += K0 * n0**2

# plt.plot(-10,0,c='k',lw=1,ls=(0,(2,1)),label=r'$K(\rho_I)$')
# plt.plot(-10,0,c='b',lw=1,ls=(0,(2,1)),label=r'$K(\rho_I)\rho^2_I$')


plt.legend(ncol=5,loc='lower center',bbox_to_anchor=[0.5,1])
plt.axis([-n.max()*0.01,n.max(),exc.min(),np.abs(exc.min())])
plt.xlabel('n')

# plt.savefig('E0.png',format='png',dpi=300)
plt.show()