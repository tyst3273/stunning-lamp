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
        self.optional_variables = ['task'] # default

        # list of mandatory variables
        self.mandatory_variables = ['prim_scale',
                                    'prim_vectors',
                                    'num_basis_atoms',
                                    'basis_atom_positions']

        # allowed variables is union of optional and mandatory variables
        self.allowed_variables = self.optional_variables+self.mandatory_variables

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

        # -----------------------------------------------------------------------------------------
        # now get the variables. if new vars are added that depend on previous variables, be careful
        # to add them in the logical order
        # -----------------------------------------------------------------------------------------

    # -----------------------------------------------------------------------------------------------

    # -----------------------------------------------------------------------------------------------
    # private methods go below here 
    # -----------------------------------------------------------------------------------------------

    def _preprocess_input(self):

        """
        strip all the blank lines, comments, and check that all mandatory variables are defined. 
        """

        duplicates = [] # store duplicate variables
        trimmed_text = [] # stripped of comments and blank lines
        for line in self.input_text:

            # if comment line or blank, skip. otherwise, process it further
            if len(line.split()) == 0 or line.strip().startswith('#'):
                continue
            else: 
     
                # strip comment off end of line
                tmp_line = line.split('#')[0].strip()

                # get the key word from the lines
                key_word = tmp_line.split('=')

                # see if it is a key word of a variable definition broken on another line
                if len(key_word) != 1:

                    # get the keyword from the list
                    key_word = key_word[0].strip()
                    
                    # see if it is an allowed keyword
                    if key_word not in self.allowed_variables:
                        message = f'key word \'{key_word}\' is not one of the allowed keywords'
                        raise_error(message)

                    # if it is mandatory, pop it out of mandatory list
                    if key_word in self.mandatory_variables:
                        index = self.mandatory_variables.index(key_word)
                        self.mandatory_variables.pop(index)

                    # check if key_word is duplicated
                    if key_word in duplicates:
                        message = f'key word \'{key_word}\' appears more than once in the input file'
                        raise_error(message)

                    # add to list to check for duplicates
                    duplicates.append(key_word)

                # append line of txt to what we want to keep
                trimmed_text.append(tmp_line)
        
        # now check that all mandatory variables are specified
        if len(self.mandatory_variables) != 0:
            message = 'the following mandatory variables are not in the input file:\n'
            for variable in self.mandatory_variables:
                message = message+f'  \'{variable}\'\n'
            raise_error(message)

    # -----------------------------------------------------------------------------------------------

    def _postprocess_input(self):

        """
        do necessary manipulations and and post-checking here
        """ 
        pass 

    # -----------------------------------------------------------------------------------------------

    # -----------------------------------------------------------------------------------------------
    # these are methods to get the variables from the input text
    # -----------------------------------------------------------------------------------------------

    # -----------------------------------------------------------------------------------------



