import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(5,10),gridspec_kw={'hspace':0.05})


base = np.loadtxt('base/constQ.txt')
base = np.nan_to_num(base,nan=0)

room = np.loadtxt('300K/constQ.txt')
room = np.nan_to_num(room,nan=0)


ax.errorbar(base[:,0], base[:,1]+0.09, yerr=base[:,2], ecolor='k', elinewidth=3, capsize=5, barsabove=True,
    markersize=8,marker='o',mfc='k',mec='k',ls=':',lw=2,color='k')
ax.errorbar(room[:,0], room[:,1], yerr=room[:,2], ecolor='r', elinewidth=3, capsize=5, barsabove=True,
    markersize=8,marker='o',mfc='r',mec='r',ls=':',lw=2,color='r')



lims = [50,100,0,0.2]
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
ax.minorticks_on()
ax.tick_params(which='both',width=1,labelsize='large')
ax.tick_params(which='major',length=5)
ax.tick_params(which='minor',length=2)
ax.axis(lims)

#ax[0].set_yticks([-10,-8,-6,-4,-2,0,2,4,6,8,10])
ax.set_ylabel('intensity (arb. units)',labelpad=5,fontsize='x-large')

ax.set_xlabel('E (meV)',labelpad=5,fontsize='x-large')

#ax[0].set_title(r'S($\bf{Q}$,E)',fontsize='x-large')

plt.savefig('constQ.png',dpi=150,bbox_inches='tight')
#plt.show()




