
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
crystal.freeze_rotation_in_supercell(coord=[0.33333333,0.33333333,0.33333333],
                                     nn_dist=5,euler_angles=[12,0,0])
crystal.write_poscar(which='supercell')
crystal.write_xyz(which='supercell')


