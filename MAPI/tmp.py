import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches


basis = np.array([[0.5,0.5],
                  [1.0,0.0],
                  [0.0,0.0],
                  [0.0,1.0],
                  [1.0,1.0]])

nat = basis.shape[0]


pos_1 = np.copy(basis)
pos_2 = np.copy(basis)

phi = -10
phi = phi*np.pi/180
R = np.array([[ np.cos(phi),-np.sin(phi)],
              [ np.sin(phi), np.cos(phi)]])

for ii in range(nat):
    pos_1[ii,:] = np.matmul(R,pos_1[ii,:])
    pos_2[ii,:] = np.matmul(R.T,pos_2[ii,:])


pos_2[:,0] = pos_2[:,0]+pos_1[-1,0]
pos_2[:,1] = pos_2[:,1]+pos_1[-1,1]


fig, ax = plt.subplots(figsize=(6,6))

ax.scatter(pos_1[:,0],pos_1[:,1],s=50,c='r')
ax.scatter(pos_2[:,0],pos_2[:,1],s=50,c='b')


r_1 = patches.Polygon(pos_1[1:,:],color='r',alpha=0.5)
r_2 = patches.Polygon(pos_2[1:,:],color='b',alpha=0.5)
ax.add_patch(r_1)
ax.add_patch(r_2)

plt.show()




