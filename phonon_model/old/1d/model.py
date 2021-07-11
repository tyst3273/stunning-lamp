import numpy as np

a = 1
pos = np.array([0,0.25])
cart = np.copy(pos)*a
num_atoms = pos.shape[0]
masses = np.array([5,2])

sc = 5
sc_a = sc*a
sc_num_atoms = sc*num_atoms
sc_pos = np.zeros(sc_num_atoms)
sc_vec = np.zeros(sc)

for ii in range(sc):
    sc_pos[ii*num_atoms:(ii+1)*num_atoms] = pos+ii
    sc_vec[ii] = ii

sc_pos = sc_pos/sc
sc_cart = sc_a*sc_pos

# -------------------------------------------------

shift = -sc*(sc_vec > sc/2).astype(int)
shift = shift+sc*(sc_vec <= -sc/2).astype(int)
sc_vec = sc_vec+shift
sc_vec = sc_vec*a

sc_dist = np.abs(sc_vec)

# -------------------------------------------------

# if there are Np atoms in primitive cell, there 
# are Np*(Np+1)/2 independent dynmat elements

num_pairs = num_atoms*(num_atoms+1)//2

# hard coded for basis size
k = np.array([1,2,1])
xi = np.array([2,1,1])

dyn_mat = np.zeros((num_atoms,num_atoms),dtype=complex)

count = 0
for ii in range(num_atoms):
    for jj in range(ii,num_atoms):

        dist = np.abs(pos[ii]-pos[jj])

        print(sc_dist+dist)
#        fc = -k[count]*np.exp(-sc_dist

        count = count+1
















