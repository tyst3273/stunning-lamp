
# -----------------------------------------------------------------------------------------------

import numpy as np

# -----------------------------------------------------------------------------------------------

import mod_lattice 
import mod_potential
import mod_io

# -----------------------------------------------------------------------------------------------

# initialize and print lattice info

p_lat_vecs = [[0,1.98918,1.98918],[1.98918,0,1.98918],[1.98918,1.98918,0]]
p_basis_pos = [[0,0,0]]
sc_matrix= [[2,0,0],[0,2,0],[0,0,2]]

lat = mod_lattice.lattice(p_lat_vecs=p_lat_vecs,p_basis_pos=p_basis_pos,sc_matrix=sc_matrix)
#print(lat.sc_rel_cart)

# -----------------------------------------------------------------------------------------------

# initialze potential 
fc = mod_potential.from_file(lat,fc_file='fc.dat')
fc.q_path(qmin=[0,0,0],qmax=[2,0,0],num_q=101)
fc.calculate(lat)

# -----------------------------------------------------------------------------------------------

# write some info to various output files
mod_io.write_lattice(lat)
mod_io.write_nn(lat)

# -----------------------------------------------------------------------------------------------

import matplotlib.pyplot as plt
for bb in range(3):
    plt.plot(fc.q_points,fc.modes[:,bb])
plt.show()





































