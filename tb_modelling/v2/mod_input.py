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

    def parse_input_file(self,input_file='input_params'):

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
        strip all the blank lines, comments, and check that all mandatroy variables are defined. 
        """

        duplicates = [] # store duplicate variables
        trimmed_text = [] # stripped of comments and blank lines
        for line in self.input_text:

            # if comment line or blank, skip. 
            if len(line.split()) == 0 or line.strip().startswith('#'):
                continue

            # otherwise, process the line further
            else: 
     
                # strip comment off end of line
                tmp_line = line.split('#')[0].strip()

                print(tmp_line)

                # see if it is an allowed keyword
                if tmp_line in self.allowed_variables:

                    # if it is mandatory, pop it out of mandatory list
                    if tmp_line in self.mandatory_variables:
                        index = self.mandatory_variables.index(tmp_line)
                        self.mandatory_variables.pop(index)

                    # check if tmp_line is duplicated
                    if tmp_line in duplicates:
                        message = f'key word \'{tmp_line}\' appears more than once in the input file'
                        raise_error(message)

                    # add to list to check for duplicates
                    duplicates.append(tmp_line)

                # append line of txt to what we want to keep
                trimmed_text.append(tmp_line)
        
        # now check that all mandatory variables are specified
        if len(self.mandatory_variables) != 0:
            message = 'the following mandatory variables are not in the input file:\n'
            for variable in self.mandatory_variables:
                message = message+f'  \'{variable}\'\n'
            raise_error(message)

        # -------------------------------------------------------------------------------------------
        # now go put all the variables and thier params into a dictionary
        # -------------------------------------------------------------------------------------------

        self.input_dict = {} # holds the input variables 
        kw_inds = [] # holds inds of key words in the trimmed_text list

        # go get the inds of the key words from the file
        for ii in range(len(trimmed_text)):
            line = trimmed_text[ii]
            if line in self.allowed_variables:
                kw_inds.append(ii)

        # now put the variables in the dictionary
        for ii in range(len(kw_inds)-1):
            ind = kw_inds[ii]

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



