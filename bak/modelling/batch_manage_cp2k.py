import os
import sys
import numpy as np

# --------------------------------------------------------------------------------------------------

allowed_jobs = ['setup','run','dry_run','dfset','none','get_energy']

job = 'none'
if len(sys.argv) != 1:
    job = sys.argv[1]
if job not in allowed_jobs:
    exit(f'job: \'{job}\' unknown')

ls = os.listdir()
dirs = []
for _ in ls:
    if _.startswith('cp2k_job'):
        dirs.append(_)

if job == 'none':
    print('allowed jobs:\n',allowed_jobs)

# --------------------------------------------------------------------------------------------------

if job == 'setup':

    with open('template.inp','r') as _:
        template = _.readlines()

    for job_dir in dirs:
        with open(os.path.join(job_dir,'cp2k.pos'),'r') as _:
            pos = _.read()
        with open(os.path.join(job_dir,'cp2k.inp'),'w') as _:
            for line in template:
                if line.strip() == '!COORDS':
                    _.write(pos)
                    _.write('\n')
                else:
                    _.write(line)

# --------------------------------------------------------------------------------------------------

if job in ['run','dry_run']:
    
    run_cmd = '/usr/bin/mpirun -np 16 cp2k.popt -i cp2k.inp -o cp2k.out'
    cwd = os.getcwd()
    print(f'cwd: {cwd}\n')

    for job_dir in dirs:
        os.chdir(job_dir)
        print(os.getcwd())
        print(f' {run_cmd}\n')

        if job == 'run':
            val = os.system(run_cmd)
            if val != 0:
                print(f'fuck! job in:\n {os.getcwd()}\n' + \
                      'failed continuining to next job with fingers crossed...\n')

        os.chdir(cwd)

# --------------------------------------------------------------------------------------------------

if job == 'dfset':

    force_conv = 1/2 # ha/bohr to ry/bohr
    len_conv = 1.88973 # angstrom to borh

    # get undistorted positions
    with open('cp2k.undistorted','r') as _:
        _.readline()
        lat_vecs = np.zeros((3,3),dtype=float)
        lat_vecs[:,0] = _.readline().strip().split()[1:]
        lat_vecs[:,1] = _.readline().strip().split()[1:]
        lat_vecs[:,2] = _.readline().strip().split()[1:]
        print('lat_vecs:\n',lat_vecs)
        for ii in range(4):
            _.readline()
        _pos = _.readlines()
        n_at = len(_pos)-1
        pos = np.zeros((n_at,3),dtype=float)
        for ii in range(n_at):
            pos[ii,:] = _pos[ii].strip().split()[1:]
            pos[ii,:] = np.matmul(lat_vecs,pos[ii,:])
    undistorted = pos*len_conv # cartesian coords in Bohr

    f_out = open('DFSET','w')

    for job_dir in dirs:

        # get displaced positions
        with open(os.path.join(job_dir,'cp2k.pos'),'r') as _:
            for ii in range(4):
                _.readline()
            for ii in range(4):
                _.readline()
            _pos = _.readlines()
            pos = np.zeros((n_at,3),dtype=float)
            for ii in range(n_at):
                pos[ii,:] = _pos[ii].strip().split()[1:]
                pos[ii,:] = np.matmul(lat_vecs,pos[ii,:])            
        pos = pos*len_conv # cartesian coords in Bohr
        pos = pos-undistorted
        
        # find force file
        files = os.listdir(job_dir)
        for ff in files:
            if ff.endswith('xyz'):
                force_file = ff
                break
    
        # get forces
        with open(os.path.join(job_dir,force_file),'r') as _:
            forces = np.zeros((n_at,3))
            for ii in range(4):
                _.readline()
            _forces = _.readlines()
            for ii in range(n_at):
                forces[ii,:] = _forces[ii].strip().split()[3:]
        forces = forces*force_conv

        # now write to DFSET file
        f_out.write(f'# {job_dir} \n')
        for ii in range(n_at):
            f_out.write(f'{pos[ii,0]: 14.12f} {pos[ii,1]: 14.12f} {pos[ii,2]: 14.12f}  '\
                        f'{forces[ii,0]: 14.12f} {forces[ii,1]: 14.12f} {forces[ii,2]: 14.12f}\n')

    f_out.close()
            
# --------------------------------------------------------------------------------------------------

if job == 'get_energy':

    for job_dir in dirs:
        ls = os.listdir(job_dir)
        log_file = None
        for ff in ls:
            if ff.endswith('.out'):
                log_file = ff
                break
        if log_file == None:
            print(f'no log file in {job_dir}. skipping!')
            continue

        with open(os.path.join(job_dir,log_file),'r') as _:
            log_data = _.readlines()

        energies = []
        for line in log_data:
            if line.strip().startswith('- Atoms:'):
                num_atoms = int(line.strip().split()[-1])
            if line.strip().startswith('ENERGY|'):
                energies.append(float(line.strip().split()[-1]))
        num_steps = len(energies)
        energy = energies[-1]
        print(f'\n{job_dir}')
        print(f'\tnum. optimizations in file: {num_steps}')
        print(f'\tfinal energy: {energy}')
        print(f'\tnum. atoms: {num_atoms}')
        print(f'\tenergy/atom: {energy/num_atoms}')

# --------------------------------------------------------------------------------------------------



