
import numpy as np
import matplotlib.pyplot as plt


def get_v(t,w0,g):
    
    mask = 1 #t > 0.001).astype(float)

    if g < w0:
        o = np.sqrt(w0**2-g**2)
        vt = v0*np.exp(-g*t)*(np.cos(o*t)-g/o*np.sin(o*t))*mask
        xt = v0*np.exp(-g*t)*np.sin(o*t)/o*mask
    elif g > w0:
        o = np.sqrt(g**2-w0**2)
        vt = v0*np.exp(-g*t)*(np.cosh(o*t)-g/o*np.sinh(o*t))*mask
        xt = v0*np.exp(-g*t)*np.sinh(o*t)/o*mask
    else:
        o = 0
        vt = v0*np.exp(-g*t)*(1-g*t)*mask
        xt = v0*np.exp(-g*t)*t*mask

    ind = np.flatnonzero(np.abs(vt) < 0.001)[0]
    print(ind)


    vrms = np.sqrt(np.trapz(vt[:ind]**2,t[:ind])/2/t[:ind].max())
    print(vrms,2*vrms/xt[:ind].max())

    return o, vt, xt, ind


v0 = 1
w0 = 5


nt = 10001
t_max = np.pi*5/w0
t = np.linspace(0.0001,t_max,nt)

for gg in np.arange(0,100,0.5):
    get_v(t,w0,gg)


g = 0.01
o, vt, xt, ind  = get_v(t,w0,g)
plt.plot(t[:ind],vt[:ind],c='r',label=r'$\gamma$='+f'{g}')
plt.plot(t[:ind],xt[:ind],c='r',ls='--')
plt.show()

g = 1
o, vt, xt, ind  = get_v(t,w0,g)
plt.plot(t[:ind],vt[:ind],c='r',label=r'$\gamma$='+f'{g}')
plt.plot(t[:ind],xt[:ind],c='r',ls='--')
plt.show()

g = 5
o, vt, xt, ind  = get_v(t,w0,g)
plt.plot(t[:ind],vt[:ind],c='r',label=r'$\gamma$='+f'{g}')
plt.plot(t[:ind],xt[:ind],c='r',ls='--')
plt.show()

g = 20
o, vt, xt, ind  = get_v(t,w0,g)
plt.plot(t[:ind],vt[:ind],c='r',label=r'$\gamma$='+f'{g}')
plt.plot(t[:ind],xt[:ind],c='r',ls='--')
plt.show()

g = 100
o, vt, xt, ind  = get_v(t,w0,g)
plt.plot(t[:ind],vt[:ind],c='r',label=r'$\gamma$='+f'{g}')
plt.plot(t[:ind],xt[:ind],c='r',ls='--')
plt.show()



plt.legend()
plt.show()


