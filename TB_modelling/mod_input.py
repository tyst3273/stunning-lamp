import os

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
        
        # list of variables with defaults
        self.default_variables = ['task','energy_cutoff']

        # list of mandatory variables
        self.mandatory_variables = ['lattice_vectors',
                                    'atom_positions',
                                    'orbital_types',
                                    '_orbital_params_'] # figure out what params to assign

        # default variables 
        self.energy_cutoff = 250 # eV, integer

    # ----------------------------------------------------------------------------------------------

    def parse_input_file(self,input_file='tb_modelling.in'):

        """
        read the input file 
        """

        # check if input file exists, if not exit
        if not os.path.exits(input_file):
            message = 'input file \'{input_file}\' not found'
            raise_error(message)

        # get the text from the input file
        with open(input_file,'r') as f_input:
            self.input_txt = f_input.readlines() # read all the lines

        # now format the lines and check that all mandatory variables exist
        self._preprocess_input(self)

    # -----------------------------------------------------------------------------------------------
    # private variables go below here 
    # -----------------------------------------------------------------------------------------------

    def _preprocess_input(self):

        """
        strip all the blank lines, comments, and check that all mandatory variables are defined
        """


        pass 

