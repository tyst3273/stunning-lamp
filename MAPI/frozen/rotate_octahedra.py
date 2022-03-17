
import classes


basis_pos = [[0.5, 0.5, 0.5],
             [0.0, 0.0, 0.0],
             [0.5, 0.0, 0.0],
             [0.0, 0.5, 0.0],
             [0.0, 0.0, 0.5]]
basis_types = ['Cs','Pb','Br','Br','Br']
lattice_vectors = [[6.02, 0.00, 0.00],
                   [0.00, 6.02, 0.00],
                   [0.00, 0.00, 6.02]]
crystal = classes.crystal(basis_pos,basis_types,lattice_vectors)


supercell_matrix = [3,3,3]
crystal.make_supercell(supercell_matrix)


