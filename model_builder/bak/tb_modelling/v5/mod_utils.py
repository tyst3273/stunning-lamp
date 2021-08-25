
# -------------------------------------------------------------------------------------------------

def print_message(message,message_type='MESSAGE'):

    """
    print stuff to the screen using a standard formatting
    """

    print(f'\n ** {message_type} **\n')
    print(f' {message}\n',flush=True)

# -------------------------------------------------------------------------------------------------

def raise_error(error_message):

    """
    print error message and exit
    """

    print_message(error_message,message_type='ERROR')
    exit()

# -------------------------------------------------------------------------------------------------


