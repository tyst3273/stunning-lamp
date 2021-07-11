import numpy as np

# -----------------------------------------------------------------------------------------------

def write_lattice(lattice):

    with open('lattice.out','w') as fout:

        # supercell transformation
        fout.write('\n-------------------------------------------------------\n')
        fout.write('\n\n*** primitive to supercell transformation (A=a*T) ***\n\n')
        message = ('! transformation matrix (T):\n'
                  f'  {lattice.sc_matrix[0,0]: 2.3f} {lattice.sc_matrix[0,1]: 2.3f}'
                  f' {lattice.sc_matrix[0,2]: 2.3f}\n  {lattice.sc_matrix[1,0]: 2.3f}'
                  f' {lattice.sc_matrix[1,1]: 2.3f} {lattice.sc_matrix[1,2]: 2.3f}\n'
                  f'  {lattice.sc_matrix[2,0]: 2.3f} {lattice.sc_matrix[2,1]: 2.3f}'
                  f' {lattice.sc_matrix[2,2]: 2.3f}\n')
        fout.write(message)
        message = ('! inverse of transformation matrix (T^-1):\n'
                  f'  {lattice.sc_matrix_inv[0,0]: 2.3f} {lattice.sc_matrix_inv[0,1]: 2.3f}'
                  f' {lattice.sc_matrix_inv[0,2]: 2.3f}\n  {lattice.sc_matrix_inv[1,0]: 2.3f}'
                  f' {lattice.sc_matrix_inv[1,1]: 2.3f} {lattice.sc_matrix_inv[1,2]: 2.3f}\n'
                  f'  {lattice.sc_matrix_inv[2,0]: 2.3f} {lattice.sc_matrix_inv[2,1]: 2.3f}'
                  f' {lattice.sc_matrix_inv[2,2]: 2.3f}\n')
        fout.write(message)

        # write primitive cell info
        fout.write('\n-------------------------------------------------------\n')
        fout.write('\n\n*** primitive cell (a) ***\n')
        fout.write(f' there are {lattice.p_num_atoms} atoms in the primitive cell\n\n')
        message = '! basis positions (lattice coords)\n'
        for atom in range(lattice.p_num_atoms):
            message = message+(f'  [{atom+1}]    type: {lattice.p_basis_types[atom]}'
                f'    pos: {lattice.p_basis_pos[atom,0]: 2.3f} {lattice.p_basis_pos[atom,1]: 2.3f}'
                f' {lattice.p_basis_pos[atom,2]: 2.3f}\n')
        fout.write(message)
        message = '! basis positions (cartesian coords, Angstrom)\n'
        for atom in range(lattice.p_num_atoms):
            message = message+(f'  [{atom+1}]    type: {lattice.p_basis_types[atom]}'
                f'    pos: {lattice.p_basis_cart[atom,0]: 2.3f} {lattice.p_basis_cart[atom,1]: 2.3f}'
                f' {lattice.p_basis_cart[atom,2]: 2.3f}\n')
        fout.write(message)
    
        message = ('\n! real space lattice (a) (Angstrom):\n'
                   f'  {lattice.p_lat_vecs[0,0]: 2.3f} {lattice.p_lat_vecs[0,1]: 2.3f}'
                   f' {lattice.p_lat_vecs[0,2]: 2.3f}\n  {lattice.p_lat_vecs[1,0]: 2.3f}'
                   f' {lattice.p_lat_vecs[1,1]: 2.3f} {lattice.p_lat_vecs[1,2]: 2.3f}\n'
                   f'  {lattice.p_lat_vecs[2,0]: 2.3f} {lattice.p_lat_vecs[2,1]: 2.3f}'
                   f' {lattice.p_lat_vecs[2,2]: 2.3f}\n')
        fout.write(message)
        message = ('! reciprocal space lattice (1/Angstrom)):\n'
                    f'  {lattice.p_recip_lat_vecs[0,0]: 2.3f} {lattice.p_recip_lat_vecs[0,1]: 2.3f}'
                    f' {lattice.p_recip_lat_vecs[0,2]: 2.3f}\n  {lattice.p_recip_lat_vecs[1,0]: 2.3f}'
                    f' {lattice.p_recip_lat_vecs[1,1]: 2.3f} {lattice.p_recip_lat_vecs[1,2]: 2.3f}\n'
                    f'  {lattice.p_recip_lat_vecs[2,0]: 2.3f} {lattice.p_recip_lat_vecs[2,1]: 2.3f}'
                    f' {lattice.p_recip_lat_vecs[2,2]: 2.3f}\n')
        fout.write(message)

        message = '\n! inverse of primitive lattice vectors (a^-1)\n'
        message = message+(f'  {lattice.p_lat_inv[0,0]: 2.3f} {lattice.p_lat_inv[0,1]: 2.3f}'
                           f' {lattice.p_lat_inv[0,2]: 2.3f}\n')
        message = message+(f'  {lattice.p_lat_inv[1,0]: 2.3f} {lattice.p_lat_inv[1,1]: 2.3f}'
                           f' {lattice.p_lat_inv[1,2]: 2.3f}\n')
        message = message+(f'  {lattice.p_lat_inv[2,0]: 2.3f} {lattice.p_lat_inv[2,1]: 2.3f}'
                           f' {lattice.p_lat_inv[2,2]: 2.3f}\n')
        fout.write(message)

        # write supercell lattice info
        fout.write('\n-------------------------------------------------------\n')
        fout.write('\n\n*** supercell (A) ***\n')
        fout.write(f' there are {lattice.sc_num_atoms} atoms in the supercell\n')
        message = '! sc basis positions (lattice coords)\n'
        for atom in range(lattice.sc_num_atoms):
            message = message+(f'  [{atom+1}]    type: {lattice.sc_basis_types[atom]}'
                f'    pos: {lattice.sc_basis_pos[atom,0]: 2.3f} {lattice.sc_basis_pos[atom,1]: 2.3f}'
                f' {lattice.sc_basis_pos[atom,2]: 2.3f}\n')
        fout.write(message)
        message = '! sc basis positions (cartesian coords, Angstrom)\n'
        for atom in range(lattice.sc_num_atoms):
            message = message+(f'  [{atom+1}]    type: {lattice.sc_basis_types[atom]}'
                f'    pos: {lattice.sc_basis_cart[atom,0]: 2.3f} {lattice.sc_basis_cart[atom,1]: 2.3f}'
                f' {lattice.sc_basis_cart[atom,2]: 2.3f}\n')
        fout.write(message)

        message = ('\n! sc real space lattice (A) (Angstrom):\n'
                    f'  {lattice.sc_lat_vecs[0,0]: 2.3f} {lattice.sc_lat_vecs[0,1]: 2.3f}'
                    f' {lattice.sc_lat_vecs[0,2]: 2.3f}\n  {lattice.sc_lat_vecs[1,0]: 2.3f}'
                    f' {lattice.sc_lat_vecs[1,1]: 2.3f} {lattice.sc_lat_vecs[1,2]: 2.3f}\n'
                    f'  {lattice.sc_lat_vecs[2,0]: 2.3f} {lattice.sc_lat_vecs[2,1]: 2.3f}'
                    f' {lattice.sc_lat_vecs[2,2]: 2.3f}\n')
        fout.write(message)
        message = ('! sc reciprocal space lattice (1/Angstrom)):\n'
                    f'  {lattice.sc_recip_lat_vecs[0,0]: 2.3f} {lattice.sc_recip_lat_vecs[0,1]: 2.3f}'
                    f' {lattice.sc_recip_lat_vecs[0,2]: 2.3f}\n  {lattice.sc_recip_lat_vecs[1,0]: 2.3f}'
                    f' {lattice.sc_recip_lat_vecs[1,1]: 2.3f} {lattice.sc_recip_lat_vecs[1,2]: 2.3f}\n'
                    f'  {lattice.sc_recip_lat_vecs[2,0]: 2.3f} {lattice.sc_recip_lat_vecs[2,1]: 2.3f}'
                    f' {lattice.sc_recip_lat_vecs[2,2]: 2.3f}\n')
        fout.write(message)

        message = '\n! inverse of supercell lattice vectors (A^-1)\n'
        message = message+(f'  {lattice.sc_lat_inv[0,0]: 2.3f} {lattice.sc_lat_inv[0,1]: 2.3f}'
                           f' {lattice.sc_lat_inv[0,2]: 2.3f}\n')
        message = message+(f'  {lattice.sc_lat_inv[1,0]: 2.3f} {lattice.sc_lat_inv[1,1]: 2.3f}'
                           f' {lattice.sc_lat_inv[1,2]: 2.3f}\n')
        message = message+(f'  {lattice.sc_lat_inv[2,0]: 2.3f} {lattice.sc_lat_inv[2,1]: 2.3f}'
                           f' {lattice.sc_lat_inv[2,2]: 2.3f}\n')
        fout.write(message)

# -----------------------------------------------------------------------------------------------

def write_nn(lattice):

    with open('coordination.out','w') as fout:
        
        fout.write('! distance in Angstrom and coordination\n')
        for atom in range(lattice.p_num_atoms):

            dist = np.round(lattice.sc_rel_sph[atom,:,0],6)
            unique, count = np.unique(dist,return_counts=True)
            message_1 = f'  [{atom+1}]:    dist:'
            message_2 = f'         coord:'
            for ii in range(1,len(unique)):
                message_1 = message_1+f' {unique[ii]:2.3f}'
                message_2 = message_2+f' {count[ii]:4g}'
            fout.write(message_1+'\n')
            fout.write(message_2+'\n')
    
# ----------------------------------------------------------------------------------------------

def write_pairs(lattice,pot):

    with open('pairs.out','w') as fout:

        fout.write('! pair types for 2-site potential matrix\n\n')
        for pair in range(pot.num_pair_types):
            atom_types = pot.pair_atoms[pair]
            fout.write(f'! pair type: {pair}; atom types: {atom_types[0]:g}, {atom_types[1]:g}\n')
            message = '  '
            for atom in range(lattice.sc_num_atoms):
                for ii in range(lattice.sc_num_atoms):
                    message = message+f' {pot.pair_matrix[pair][atom,ii]:g}'
                message = message+'\n  '
            fout.write(message+'\n')

# ----------------------------------------------------------------------------------------------















