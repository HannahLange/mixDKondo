import numpy as np
## Helper functions


def log(string):
    """
    Print progress
    """

    print("-"*20 + " " + string + " ..." + "-"*20)



def xy_to_snake(x,y,width):
    """
    Maps [x,y] lattice coordinates to 1D snake
    """

    return x*width + y
