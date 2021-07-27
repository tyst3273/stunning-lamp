"""
various utility stuff like wrappers to print to screen, raise error, etc
"""

# -----------------------------------------------------------------------------------------------------

def print_stdout(message,msg_type='NOTE'):

    """
    print to the screen using common formatting
    """
    
    print(f'\n ** {msg_type} **')
    print(f' {message}\n',flush=True)

# -----------------------------------------------------------------------------------------------------

def raise_error(message):

    print_stdout(message,msg_type='ERROR')
    exit()

# -----------------------------------------------------------------------------------------------------


