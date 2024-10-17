import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl

mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = r"\usepackage{bm}"

nx = 501
ny = 501
x, y = np.meshgrid(np.linspace(-1,1,nx),np.linspace(-1,1,ny),indexing='ij')

fig, ax = plt.subplots(2,3,figsize=(8,5.25),gridspec_kw={'wspace':0.075,'hspace':0.175})

#tau_x = np.array([0.0,0.5,0.0])
#tau_y = np.array([0.0,0.0,0.5])
#s = [100,50,50]
#c = ['k','w','w']
#b = [7,5,5]

tau_x = np.array([0.0,0.5])
tau_y = np.array([0.0,0.5])
s = [100,100]
c = ['k','w']
b = [1,1]

extent = [-1,1,-1,1]

Gx = 1
Gy = 0
exp_iQr = np.real(np.exp(2j*np.pi*(Gx*x+Gy*y)))
ax[0,0].imshow(exp_iQr.T,cmap='coolwarm',aspect='auto',
                interpolation='none',origin='lower',extent=extent)
ax[0,0].set_title(r'$\bm{G}$'+f'=({Gx},{Gy})',fontsize='large',y=0.975)
Fq = np.abs(np.sum(b*np.exp(2j*np.pi*(Gx*tau_x+Gy*tau_y))))
print(Fq)

Gx = 2
Gy = 0
exp_iQr = np.real(np.exp(2j*np.pi*(Gx*x+Gy*y)))
ax[0,1].imshow(exp_iQr.T,cmap='coolwarm',aspect='auto',
                interpolation='none',origin='lower',extent=extent)
ax[0,1].set_title(r'$\bm{G}$'+f'=({Gx},{Gy})',fontsize='large',y=0.975)
Fq = np.abs(np.sum(b*np.exp(2j*np.pi*(Gx*tau_x+Gy*tau_y))))
print(Fq)

Gx = 3
Gy = 0
exp_iQr = np.real(np.exp(2j*np.pi*(Gx*x+Gy*y)))
ax[0,2].imshow(exp_iQr.T,cmap='coolwarm',aspect='auto',
                interpolation='none',origin='lower',extent=extent)
ax[0,2].set_title(r'$\bm{G}$'+f'=({Gx},{Gy})',fontsize='large',y=0.975)
Fq = np.abs(np.sum(b*np.exp(2j*np.pi*(Gx*tau_x+Gy*tau_y))))
print(Fq)

Gx = 1
Gy = 1
exp_iQr = np.real(np.exp(2j*np.pi*(Gx*x+Gy*y)))
ax[1,0].imshow(exp_iQr.T,cmap='coolwarm',aspect='auto',
                interpolation='none',origin='lower',extent=extent)
ax[1,0].set_title(r'$\bm{G}$'+f'=({Gx},{Gy})',fontsize='large',y=0.975)
Fq = np.abs(np.sum(b*np.exp(2j*np.pi*(Gx*tau_x+Gy*tau_y))))
print(Fq)

Gx = 2
Gy = 2
exp_iQr = np.real(np.exp(2j*np.pi*(Gx*x+Gy*y)))
ax[1,1].imshow(exp_iQr.T,cmap='coolwarm',aspect='auto',
                interpolation='none',origin='lower',extent=extent)
ax[1,1].set_title(r'$\bm{G}$'+f'=({Gx},{Gy})',fontsize='large',y=0.975)
Fq = np.abs(np.sum(b*np.exp(2j*np.pi*(Gx*tau_x+Gy*tau_y))))
print(Fq)

Gx = 2
Gy = 1
exp_iQr = np.real(np.exp(2j*np.pi*(Gx*x+Gy*y)))
ax[1,2].imshow(exp_iQr.T,cmap='coolwarm',aspect='auto',
                interpolation='none',origin='lower',extent=extent)
ax[1,2].set_title(r'$\bm{G}$'+f'=({Gx},{Gy})',fontsize='large',y=0.975)
Fq = np.abs(np.sum(b*np.exp(2j*np.pi*(Gx*tau_x+Gy*tau_y))))
print(Fq)


for ii in range(2):
    for jj in range(3):
        
        ax[ii,jj].plot([-1,1],[0,0],lw=1,c=(0.25,0.25,0.25),ls=(0,(2,1)))
        ax[ii,jj].plot([0,0],[-1,1],lw=1,c=(0.25,0.25,0.25),ls=(0,(2,1)))

        for xx in range(-1,1):
            for yy in range(-1,1):
                ax[ii,jj].scatter(tau_x+xx,tau_y+yy,color=c,s=s,clip_on=False,edgecolors='k',
                            linewidths=1.5,zorder=1000)



for ii in range(2):
    for jj in range(3):
        for axis in ['top','bottom','left','right']:
            ax[ii,jj].spines[axis].set_linewidth(1.5)

        ax[ii,jj].minorticks_on()
        ax[ii,jj].tick_params(which='both',width=1,labelsize='x-large')
        ax[ii,jj].tick_params(which='major',length=5)
        ax[ii,jj].tick_params(which='minor',length=2)

        ax[ii,jj].axis([-1,1,-1,1])
        ax[ii,jj].set_xticks([-1,0,1])
        ax[ii,jj].set_yticks([-1,0,1])


#ax.set_xticks(qpts_verts)
#ax.set_xticklabels([r'X/M$^*$',r'$\Gamma$',r'M/$\Gamma^*$'])

#ax.set_xlabel(r'$\bm{q}$',fontsize='x-large')
ax[0,1].set_yticklabels([])
ax[0,2].set_yticklabels([])
ax[1,1].set_yticklabels([])
ax[1,2].set_yticklabels([])

for ii in range(3):
    ax[0,ii].set_xticklabels([])
    ax[1,ii].set_xlabel(r'$\tau_x$',fontsize='x-large',labelpad=0)

ax[0,0].set_ylabel(r'$\tau_y$',fontsize='x-large',rotation='horizontal',labelpad=5)
ax[1,0].set_ylabel(r'$\tau_y$',fontsize='x-large',rotation='horizontal',labelpad=5)

#ax[1].annotate(r'$\frac{k_b T}{E_g}$='+f'{T/(2*E0):3.2f}',xy=(0.35,0.5),xycoords='axes fraction',
#                fontsize='x-large')

plt.savefig('bragg.pdf',dpi=300,bbox_inches='tight')
#plt.show()



