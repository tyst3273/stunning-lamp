import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

cutoff = 0.2

x = np.linspace(-5,5,501)
y = np.linspace(-5,5,501)
print(x)
x, y = np.meshgrid(x,y,indexing='xy')

r = np.sqrt(x**2+y**2) # len of coordinate vector

# mask the r == 0 zone 
mask = (r <= cutoff).astype(float)
r_masked = r+mask

rx = x/r_masked
ry = y/r_masked

q = 1
p = np.array([0,1.0])*q # dipole vector
r_dot_p = rx*p[0]+ry*p[1] 
k_c = 1.0 # coulombs constant
Ex = (3*r_dot_p*rx-p[0])/r_masked**3
Ey = (3*r_dot_p*ry-p[1])/r_masked**3

# now set the r==0 zone to 0 to hide it
mask = (r >= cutoff).astype(float)
Ex *= mask
Ey *= mask
E = np.sqrt(Ex**2+Ey**2)

fig, ax = plt.subplots(figsize=(6,6))
norm = Normalize(vmin=-1,vmax=0.5)
plt.streamplot(x,y,Ex,Ey,broken_streamlines=False,density=0.75,color=E,cmap='Blues',norm=norm)
plt.annotate('',xy=[0,0.5],xytext=[0,-0.5],xycoords='data',
            arrowprops={'width':2,'headwidth':10,'headlength':10,'color':'k'})
plt.plot(0,0,marker='o',ms=50,c='m',mec='k',mew=2)

for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
ax.minorticks_on()
ax.tick_params(which='both',width=1,labelsize='x-large')
ax.tick_params(which='major',length=5)
ax.tick_params(which='minor',length=2)
ax.set_rasterized = True

ax.set_xticks(np.arange(-5,6))
ax.set_yticks(np.arange(-5,6))

ax.set_xlabel('x',fontsize='x-large')
ax.set_ylabel('y',fontsize='x-large')

plt.savefig('dipole.png',dpi=300,bbox_inches='tight')
plt.show()






