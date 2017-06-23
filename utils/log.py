import numpy

def log(x, b=numpy.exp(1)):
    """
    Calculates the logarithm of base b

    INPUTS
    ------
    x : array-like
    b : float, int, default e
    """
    return numpy.log(x) / numpy.log(b)
