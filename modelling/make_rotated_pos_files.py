
import modelling_tools as mtools


io = mtools.interface_io()



atom_pos = [[0.5, 0.5, 0.5],
            [0.0, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            [0.0, 0.5, 0.0],
            [0.0, 0.0, 0.5]]
atom_types = ['Cs','Pb','Br','Br','Br']
lattice_vectors = [[6.02, 0.00, 0.00],
                   [0.00, 6.02, 0.00],
                   [0.00, 0.00, 6.02]]

crystal = mtools.crystal(atom_pos,atom_types,lattice_vectors)

supercell_matrix = [3,3,3]
angles = [8,16,24]

for angle in angles:
    f_name = f'cp2k_job.{angle:g}_degrees'
    crystal.make_supercell(supercell_matrix)
    crystal.freeze_rotation_in_supercell(coord=[0.33333333,0.33333333,0.33333333],
                                         nn_dist=4,euler_angles=[angle,0,0])
    io.write_xyz(crystal,out_file=f_name+'/pos.xyz',which='supercell')
    io.write_cp2k(crystal,out_file=f_name+'/cp2k.pos',which='supercell')
    io.write_poscar(crystal,out_file=f_name+'/POSCAR',which='supercell')
    io.write_abivars(crystal,out_file=f_name+'/abivars',which='supercell')


