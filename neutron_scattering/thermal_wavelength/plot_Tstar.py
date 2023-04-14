
import numpy as np
import matplotlib.pyplot as plt


hbar = 1.054572e-34            
kb = 1.380649e-23
amu2kg = 1.660540e-27


def T_for_thermal_wavelenght(mass,lam=1e-10):
    mass = mass*amu2kg
    return 2*np.pi*hbar**2/(mass*kb*lam**2)
    
#emass = 0.00054858
#print(T_for_thermal_wavelenght(emass))

masses = np.loadtxt('masses.txt',dtype=object)
masses = masses[:,[1,3]]
natoms = masses.shape[0]

T0 = np.zeros(natoms,dtype=float)

for ii in range(natoms):
    atom = masses[ii,0]
    mass = float(masses[ii,1])
    T0[ii] = T_for_thermal_wavelenght(mass)


fig, ax = plt.subplots(figsize=(9,6))

inds = np.arange(natoms)
ax.plot(inds,T0,marker='o',ms=10,c='k',ls=':',lw=1)

print(T0)

"""
for ii in range(natoms):
    T = T0[ii]
    if ii % 2 == 0 or ii < 8:
        yshift = 0.1*T
        xshift = 0.1
    else:
        xshift = -0.1
        yshift = -0.25*T
    ax.annotate(masses[ii,0],color='k',xycoords='data',xy=(ii+xshift,T+yshift),
                   fontsize='x-large',fontweight='bold')
"""

for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
ax.minorticks_on()
ax.tick_params(axis='x', which='minor', bottom=False)
ax.tick_params(which='both',width=1,labelsize='large')
ax.tick_params(which='major',length=5)
ax.tick_params(which='minor',length=2)

ax.set_xticks(inds)
ax.set_xticklabels(masses[:,0])
ax.set_ylabel(r'T$^*$ [K]',fontsize='x-large',rotation='vertical',labelpad=5)
ax.set_yscale('log')

ax.axis([-1,23,5,350])

plt.savefig('Tstar.pdf',bbox_inches='tight')

plt.show()
    




