import os

def mkdirp(path):
    """
    This function makes a directory if it doesn't exist, otherwise it doesn't do anything.
    It is a python wrapper for the "mkdir -p" command in bash
    """
    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)
