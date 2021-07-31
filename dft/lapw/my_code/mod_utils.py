
"""
various utility stuff like wrappers to print to screen, raise errors, etc
"""

# -----------------------------------------------------------------------------------------------------

def print_message(message,message_type='NOTE'):

    """
    print to the screen using common formatting
    """
    
    print(f'\n ** {message_type} **')
    print(f' {message}\n',flush=True)

# -----------------------------------------------------------------------------------------------------

def raise_error(message):

    print_message(message,message_type='ERROR')
    exit()

# -----------------------------------------------------------------------------------------------------


