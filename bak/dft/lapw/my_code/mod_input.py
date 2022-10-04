# system modules
import numpy as np
import os

# custom modules
from mod_utils import print_message, raise_error

# ----------------------------------------------------------------------------------------------------

class input_variables:

    """
    initialize default varibales, parse the file, and check that variables make sense
    """

    # -------------------------------------------------------------------------------------------------

    def __init__(self):

        """
        set the defaults and hold all the key words
        """

        # key words for common variables
        self.common_variables = ['tasks',
                                 'atom_positions',
                                 'atom_types',
                                 'lattice_vectors',
                                 'lattice_scale',
                                 'kpoints',
                                 'kpoints_type']

        # key words for dft variables
        self.dft_variables = ['energy_cutoff']

        # key words for tight binding variables
        self.tb_variables = []

        # to check that all given key words are allowed 
        self.allowed_variables = self.common_variables+self.dft_variables+self.tb_variables

        # set the default variables
        self.tasks = ['scf'] # what task to do
        self.energy_cutoff = 250 # pw cutoff in eV
        self.lattice_scale = 1 # scaling factor for the lattice
        self.kpoints_type = 'reduced_grid' # how to interpret the kpts arguments

    # -------------------------------------------------------------------------------------------------

    def parse_input_file(self,input_file):

        """
        preprocess and read the input file
        """
    
        # check if the input file exists
        self.input_file = input_file
        if not os.path.exists(input_file):
            message = f'input file \'{input_file}\' not found'
            raise_error(message)
   
        # remove blank and comment lines and make sure all key words look good
        self._preprocess_file()
    
        # --------------------------------------------------------------------------------------------
        # get the variables from the dictionary and check that they make sense. if new variables are 
        # added that depend on the others, add them in the logical order. if mandatory variables are 
        # missing, code should exit with message. iterable variables should have (or be given) a 
        # value for each task in tasks. global are common to all tasks.
        # -------------------------------------------------------------------------------------------
        
        # tasks. default, iterable
        if 'tasks' in self.var_dict:
            allowed_tasks = ['scf','scf_restart','tb']
            self.tasks = self.var_dict['tasks']
            self.num_tasks = len(self.tasks)
            if self.num_tasks < 1: # make sure there is atleast one
                message = 'you must give atleast 1 task to do'
                raise_error(message)
            for task in self.tasks: # make sure theyre all allowed
                if task not in allowed_tasks:
                    message = f'task \'{task}\' is not one of the allowed tasks'
                    raise_error(message)

        # energy_cutoff. default, iterable
        if 'energy_cutoff' in self.var_dict:
            self.energy_cutoff = self.var_dict['energy_cutoff']
            if len(self.energy_cutoff) == 1:  # if only one, replicate to equal number of tasks
                self.energy_cutoff = self.energy_cutoff*self.num_tasks
            elif len(self.energy_cutoff) != self.num_tasks: # check if not == 1 that its the right number
                message = 'number of energy cutoffs should be 1 or equal to the number of tasks'
                raise_error(message)
            try: # now convert to floats
                self.energy_cutoff = [float(x) for x in self.energy_cutoff]
            except: 
                message = 'energy cutoffs should be floats'
                raise_error(message)

        # lattice_scale. default, iterable
        if 'lattice_scale' in self.var_dict:
            self.lattice_scale = self.var_dict['lattice_scale']
            if len(self.lattice_scale) == 1: # if only one, replicate to equal number of tasks
                self.lattice_scale = self.lattice_scale*self.num_tasks
            elif len(self.lattice_scale) != self.num_tasks: # check if not == 1 that its the right number
                message = 'number of lattice scales should be 1 or equal to the number of tasks'
                raise_error(message)
            try: # now convert to floats
                self.lattice_scale = [float(x) for x in self.lattice_scale]
            except:
                message = 'lattice scales should be floats'
                raise_error(message)

        # atom_positions. mandatory, global
        if 'atom_positions' in self.var_dict:
            atom_pos_list = self.var_dict['atom_positions']
            self.num_atoms = len(atom_pos_list) # number of atoms
            self.atom_positions = np.zeros((self.num_atoms,3),dtype=float)
            for ii in range(self.num_atoms): # put atoms into np.array
                try:
                    self.atom_positions[ii,:] = atom_pos_list[ii].split()[:] # make sure its 3d and float
                except:
                    message = 'atom positions seem wrong'
                    raise_error(message)
        else:
            message = 'key word \'atom_positions\' is mandatory'
            raise_error(message)

        # atom_types. mandatory, global
        if 'atom_types' in self.var_dict:
            self.atom_types = self.var_dict['atom_types']
            if len(self.atom_types) != self.num_atoms: # make sure same number as numer of positions
                message = 'number of atom types does not equal number of atom positions'
                raise_error(message)
        else:
            message = 'key word \'atom_types\' is mandatory'
            raise_error(message)

        # lattice_vectors. mandatory, global
        if 'lattice_vectors' in self.var_dict:
            lat_vecs_list = self.var_dict['lattice_vectors']
            if len(lat_vecs_list) != 3: # should be 3D row vectors
                message = 'lattice vectors should be a 3x3 matrix of floats (row vectors)'
                raise_error(message)
            self.lattice_vectors = np.zeros((3,3))
            for ii in range(3): # put into np.array
                try:
                    self.lattice_vectors[ii,:] = lat_vecs_list[ii].split()[:] # should be 3D and float
                except:
                    message = 'lattice vectors should be a 3x3 matrix of floats (row vectors)'
                    raise_error(message)
        else:
            message = 'key word \'lattice_vectors\' is mandatory'
            raise_error(message)

        # kpoints_type. default, iterable
        allowed_kpts_types = ['full_grid','reduced_grid','band_path']
        if 'kpoints_type' in self.var_dict:
            self.kpoints_type = self.var_dict['kpoints_type']
            if len(self.kpoints_type) == 1:  # if only one, replicate to equal number of tasks
                self.kpoints_type = self.kpoints_type*self.num_tasks
            elif len(self.kpoints_type) != self.num_tasks: # check if not == 1 that its the right number
                message = 'number of kpoints types should be 1 or equal to the number of tasks'
                raise_error(message)
            for kpts_type in self.kpoints_type:
                if kpts_type not in allowed_kpts_types:
                    message = f'kpoints types \'{kpts_type}\' is not one of the allowed types'
                    raise_error(message)
        
        # kpoints. mandatory, iterable
        if 'kpoints' in self.var_dict:
            self.kpoints = self.var_dict['kpoints']
            if len(self.kpoints) == 1:  # if only one, replicate to equal number of tasks
                self.kpoints = self.kpoints*self.num_tasks
            elif len(self.kpoints) != self.num_tasks: # check if not == 1 that its the right number
                message = 'number of kpoints args should be 1 or equal to the number of tasks'
                raise_error(message)
            message = 'nothing is done with the kpoints in mod_input yet!!'
            print_message(message,message_type='WARNING')

    # -------------------------------------------------------------------------------------------------
    # private methods below here
    # -------------------------------------------------------------------------------------------------

    def _preprocess_file(self):

        """
        remove blank and comment lines and strip comments, set up dictionary of key words/args
        """

        # get the text from the input file
        with open(self.input_file,'r') as f_in:
            input_text = f_in.readlines()

        # strip blank and comment lines from input text
        kw_inds = []
        trimmed_text = []
        for line in input_text:
            line = line.strip() # remove newline characters and leading spaces

            # skip blank and comment lines
            if not line or line.startswith('#'): 
                continue
            
            # remove trailing comments
            line = line.split('#')[0].strip()

            # add the trimmed line to trimmed text
            trimmed_text.append(line)
    
        # now find duplicates and kw indices
        duplicates = []
        kw_inds = []
        for ii in range(len(trimmed_text)):
            line = trimmed_text[ii]
            
            # check if its a keyword
            if line.endswith(':'): 
                kw_inds.append(ii)
                kw = line.strip(':')

                # check if its an allowed variable
                if kw not in self.allowed_variables:
                    message = f'key word \'{kw}\' is not one of the allowed variables'
                    raise_error(message)

                # check if its duplicated
                if kw in duplicates:
                    message = f'the key word \'{kw}\' appears more than once in the input file'
                    raise_error(message)
                duplicates.append(kw)

        # now put the variables in the variable dictionary
        var_dict = {}
        num_kw = len(kw_inds)
        for ii in range(num_kw-1):
            kw = trimmed_text[kw_inds[ii]].strip(':')
            var_dict[kw] = trimmed_text[kw_inds[ii]+1:kw_inds[ii+1]]
        kw = trimmed_text[kw_inds[-1]].strip(':')
        var_dict[kw] = trimmed_text[kw_inds[-1]+1:]

        self.var_dict = var_dict # dictionary containing all keywords and thier args
        
    # ----------------------------------------------------------------------------------------------------













