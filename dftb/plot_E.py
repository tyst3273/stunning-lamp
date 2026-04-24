
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

x = np.linspace(-2,2,1000)
n = np.zeros_like(x)
a = 10

n0 = (a/np.pi)**(3/2) * np.exp(-a * x**2 )

for x0 in [-2,-1,0,1,2]:
    n += (a/np.pi)**(3/2) * np.exp(-a * (x-x0)**2 )

r = rs(n)
K = ex(r)+ec(r)

# plt.plot(x,n,c='r')
# plt.plot(x,n0,c='b')

plt.plot(x,K,c='k',lw=1,label=r'$K(n_0)$')
plt.plot(x,K*n**2,c='m',lw=1,label=r'$K(n_0)n^2_0$')

plt.axhline(0,lw=0.5,ls=(0,(2,1)),c='k')

I = np.zeros_like(K)

for x0 in [-2,-1,0,1,2]:
    n0 = (a/np.pi)**(3/2) * np.exp(-a * (x-x0)**2 )
    r0 = rs(n0)
    K0 = ex(r0)+ec(r0)
    plt.plot(x,K0,c='k',lw=1,ls=(0,(2,1)))
    plt.plot(x,K0*n0**2,c='b',lw=1,ls=(0,(2,1)))
    I += K0 * n0**2

plt.plot(-10,0,c='k',lw=1,ls=(0,(2,1)),label=r'$K(\rho_I)$')
plt.plot(-10,0,c='b',lw=1,ls=(0,(2,1)),label=r'$K(\rho_I)\rho^2_I$')

plt.plot(x,I,c='m',lw=1,ls=(0,(2,1)),label=r'$\sum_I K(\rho_I)\rho^2_I$')

print(np.trapz(y=I,dx=x[1]-x[0]))
print(np.trapz(y=K*n**2,dx=x[1]-x[0]))

plt.legend(ncol=5,loc='lower center',bbox_to_anchor=[0.5,1])
plt.axis([-2,2,-1.5,0.1])
plt.xlabel('|r|')

plt.savefig('E0.png',format='png',dpi=300)
plt.show()