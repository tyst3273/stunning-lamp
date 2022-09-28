import numpy as np
import matplotlib.pyplot as plt

#cmap = 'magma'
cmap = plt.cm.magma.with_extremes(under='k')
interp = None #'nearest'
vmin = 0
vmax = 9e-3


fig, ax = plt.subplots(2,1,figsize=(5,10),gridspec_kw={'hspace':0.05})

base = np.loadtxt('base/co_base_elastic.txt')[:,0]
#base = np.loadtxt('base/co_base.txt')[:,0]
base = np.nan_to_num(base,nan=-1e6)
base = base.reshape(360,400).T

room = np.loadtxt('300K/co_300K_elastic.txt')[:,0]
#room = np.loadtxt('300K/co_300K.txt')[:,0]
room = np.nan_to_num(room,nan=-1e6)
room = room.reshape(360,400).T

im = ax[0].imshow(base,aspect='auto',cmap=cmap,origin='lower',
    interpolation=interp,extent=[-6,12,-10,10],vmin=vmin,vmax=vmax)
im = ax[1].imshow(room,aspect='auto',cmap=cmap,origin='lower',
    interpolation=interp,extent=[-6,12,-10,10],vmin=vmin,vmax=vmax)

#cbar = fig.colorbar(im,ax=[ax[0],ax[1]],extend='both',location='top')

lims = [-6,12,-10,10]
for jj in range(2):
    for axis in ['top','bottom','left','right']:
        ax[jj].spines[axis].set_linewidth(1.5)
    ax[jj].minorticks_on()
    ax[jj].tick_params(which='both',width=1,labelsize='large')
    ax[jj].tick_params(which='major',length=5)
    ax[jj].tick_params(which='minor',length=2)
    ax[jj].axis(lims)

ax[0].set_yticks([-10,-8,-6,-4,-2,0,2,4,6,8,10])
ax[0].set_ylabel('K (r.l.u.)',labelpad=-1,fontsize='x-large')

ax[0].set_yticks([-10,-8,-6,-4,-2,0,2,4,6,8,10])
ax[1].set_ylabel('K (r.l.u.)',labelpad=-1,fontsize='x-large')

ax[1].set_xticks([-6,-4,-2,0,2,4,6,8,10,12])
ax[1].set_xlabel('H (r.l.u.)',labelpad=5,fontsize='x-large')

ax[1].set_xticks([-6,-4,-2,0,2,4,6,8,10,12])
ax[0].set_xticklabels([])

ax[0].set_title(r'$\Delta$E=0 '+'(meV),    L=1 (r.l.u.)',fontsize='x-large')

plt.savefig('CO.png',dpi=150,bbox_inches='tight')
#plt.show()




