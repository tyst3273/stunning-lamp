# import system modules
import numpy as np
import sys

# import custom modules
import mod_input

# get input file from command line or use default
if len(sys.argv) == 1:
    input_file = 'tb_modelling.in'
else:
    input_file = sys.argv[1]

# read the input file
invars = mod_input.invars()
invars.parse_input_file(input_file)


