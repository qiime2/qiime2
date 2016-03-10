import os

def get_install_path():
    return os.path.abspath(os.path.split(__file__)[0])
