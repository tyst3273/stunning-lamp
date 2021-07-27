import os
from mod_utils import raise_error

class invars:

    # ------------------------------------------------------------------------------------------------

    """
    input variables for the calculation are stored here. defaults are set, input file is parsed and 
    some checking on varibles is done. 
    """

    def __init__(self):

        """
        set defaults and store some other stuff. initialize all variables here
        """
        
        # list of variables with defaults or that arent required
        self.optional_variables = ['task', # default
                                   '_symmetric_to'] # optional

        # list of mandatory variables
        self.mandatory_variables = ['lattice_scale',
                                    'lattice_vectors',
                                    'number_of_atoms',
                                    'atom_positions',
                                    'number_of_orbitals',
                                    'orbital_atoms',
                                    'orbital_labels',
                                    'orbital_params',
                                    '_orbital_i',
                                    '_orbital_j',
                                    '_symmetric_to']

        # default variables go below here
        self.task = 'solve_k_grid' # what the code is gonna do

    # ----------------------------------------------------------------------------------------------

    def parse_input_file(self,input_file='tb_modelling.in'):

        """
        read the input file 
        """

        # check if input file exists, if not exit
        if not os.path.exists(input_file):
            message = f'input file \'{input_file}\' not found'
            raise_error(message)

        # get the text from the input file
        with open(input_file,'r') as f_input:
            self.input_text = f_input.readlines() # read all the lines

        # now format the lines and check that all mandatory variables exist
        self._preprocess_input()

        # now get the variables. if new vars are added, be careful to add them in the logical order

        task = self._parse_str('task')

    # -----------------------------------------------------------------------------------------------
    # private methods go below here 
    # -----------------------------------------------------------------------------------------------

    def _preprocess_input(self):

        """
        strip all the blank lines, comments, and check that all mandatory variables are defined
        """

        duplicates = [] # store duplicate variables
        trimmed_text = [] # stripped of comments and blank lines
        for line in self.input_text:

            # if comment line or blank
            if len(line.split()) == 0 or line.strip().startswith('#'):
                continue
            else:
                
                # strip comment off end of line
                tmp_line = line.split('#')[0].strip()

                # get the key word from the lines
                key_word = tmp_line.split('=')[0].split(':')[0].strip()

                # check if key_word is duplicated
                if key_word in duplicates:
                    message = f'key word \'{key_word}\' appears more than once in the input file'
                    raise_error(message)

                # append to what we want to keep
                trimmed_text.append(tmp_line)

        # overwrite the input txt and with trimmed data
        self.input_text = trimmed_text

    # -----------------------------------------------------------------------------------------------

    def _postprocess_input(self):

        """
        do necessary manipulations and and post-checking here
        """ 
        pass

    # -----------------------------------------------------------------------------------------------
    # these are methods to get the variables from the input text
    # -----------------------------------------------------------------------------------------------

    def _parse_str(self,key_word):

        """
        get str varaible from file
        """

        return_value = None
        for line in self.input_txt:
            if line.split('=')[0].strip() == key_word:
                return_value = line.split('=')[-1]
                return_value = return_value.split('#')[0].strip()
                return_value = str(return_value)
        return return_value

    # -----------------------------------------------------------------------------------------

    def _parse_float(self,key_word):

        """
        get float varible from file
        """

        return_value = None
        for line in self.input_txt:
            if line.split('=')[0].strip() == key_word:
                return_value = line.split('=')[-1]
                return_value = return_value.split('#')[0].strip()
                try:
                    return_value = float(return_value)
                except:
                    message = f'key word \'{key_word}\' seems wrongs.'
                    raise PSF_exception(message)
        return return_value

    # -----------------------------------------------------------------------------------------

    def _parse_int(self,key_word):

        """
        get int variable from file
        """

        return_value = None
        for line in self.input_txt:
            if line.split('=')[0].strip() == key_word:
                return_value = line.split('=')[-1]
                return_value = return_value.split('#')[0].strip()
                try:
                    return_value = int(return_value)
                except:
                    message = f'key word \'{key_word}\' seems wrongs.'
                    raise PSF_exception(message)
        return return_value

    # -----------------------------------------------------------------------------------------

    def _parse_bool(self,key_word):

        """
        get bool variable from file
        """

        return_value = None
        for line in self.input_txt:
            if line.split('=')[0].strip() == key_word:
                return_value = line.split('=')[-1]
                return_value = return_value.split('#')[0].strip()
                try:
                    return_value = bool(int(return_value))
                except:
                    message = f'key word \'{key_word}\' seems wrongs.'
                    raise PSF_exception(message)
        return return_value

    # -----------------------------------------------------------------------------------------

    def _parse_int_list(self,key_word):

        """
        get list of ints from file
        """

        return_value = None
        for line in self.input_txt:
            if line.split('=')[0].strip() == key_word:
                return_value = line.split('=')[-1]
                return_value = return_value.split('#')[0].strip()
                return_value = return_value.split()
                try:
                    return_value = [int(x) for x in return_value]
                except:
                    message = f'key word \'{key_word}\' seems wrongs.'
                    raise PSF_exception(message)
        return return_value

    # -----------------------------------------------------------------------------------------

    def _parse_float_list(self,key_word):

        """
        get list of floats from file
        """

        return_value = None
        for line in self.input_txt:
            if line.split('=')[0].strip() == key_word:
                return_value = line.split('=')[-1]
                return_value = return_value.split('#')[0].strip()
                return_value = return_value.split()
                try:
                    return_value = [float(x) for x in return_value]
                except:
                    message = f'key word \'{key_word}\' seems wrongs.'
                    raise PSF_exception(message)
        return return_value

    # ------------------------------------------------------------------------------------------

    def _parse_str_list(self,key_word):

        """
        get list of ints from file
        """

        return_value = None
        for line in self.input_txt:
            if line.split('=')[0].strip() == key_word:
                return_value = line.split('=')[-1]
                return_value = return_value.split('#')[0].strip()
                return_value = return_value.split()
                try:
                    return_value = [str(x) for x in return_value]
                except:
                    message = f'key word \'{key_word}\' seems wrongs.'
                    raise PSF_exception(message)
        return return_value

    # -----------------------------------------------------------------------------------------







