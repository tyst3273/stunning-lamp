import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = r"\usepackage{bm}"


u_q = 0.25

N = 8
x = np.arange(0,N)
xf = np.linspace(0,N,101)

q = 2*np.pi*0

cos = np.cos(x)


fig, ax = plt.subplots(figsize=(8,3))

n = [0,1,4]

for ii in range(N):
    plt.plot([ii,ii],[0-0.5,len(n)-0.5],lw=1,ls='--',c=(0.2,0.2,0.2))

for ii, nn in enumerate(n):

    q = 2*np.pi*nn/N

    cos = np.cos(q*x)

    u = u_q*cos

    _s = np.ones(N)*ii

#    plt.plot(x,_s,marker='o',ms=6,mfc='w',mec='k',mew=1,lw=0)
    plt.plot(x+u,_s,marker='o',ms=8,c='b',mew=0,lw=0)

    _cos = np.cos(q*xf)/5

    plt.plot(xf,_cos+ii,marker='o',ms=0,c='k',mew=0,lw=1,ls='-')

    plt.quiver(x,_s,u,np.zeros(N),angles='xy', scale_units='xy', scale=1,units='xy',
        headwidth=2,headlength=3,headaxislength=3)

    plt.annotate(rf'q=$\frac{{2\pi}}{{a}} \frac{{{nn}}}{{{N}}}$',xycoords='data',
            xy=(7.5,ii-0.2),fontsize='large')

#plt.anotate(


for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)

    ax.minorticks_on()
    ax.tick_params(which='both',width=1,labelsize='x-large')
    ax.tick_params(which='major',length=5)
    ax.tick_params(which='minor',length=2)


ax.set_xlabel('i [unitcell label]',fontsize='x-large')
ax.set_yticklabels([])


plt.savefig('1d_harmonic_disp.png',dpi=300,bbox_inches='tight')

plt.show()

