import numpy as np
from scipy.stats import norm

def pearson(x, y, CI=None):
    """
    Calculate the pearson correlation coefficient rho, defined as

        sum((x - mu_x)*(y - mu_y)) / \
           [ sqrt(sum(x - mu_x)^2) * sqrt(sum(y - mu_y)^2)]

    Confidence intervals are also calculated using a Fisher transformation.
    https://en.wikipedia.org/wiki/Fisher_transformation

    INPUTS
    ------
    x, y : NumPy arrays
        x and y must have same length
    CI : float, default None
        confidence interval, e.g. 95%

    RETURNS
    -------
    (rho, lower, upper) : (float, float, float)
        If CI is supplied. Otherwise returns rho.
    """

    n = len(x)
    rho = np.sum( (x - np.mean(x)) * (y - np.mean(y)) ) /\
            ( np.sqrt(np.sum((x-np.mean(x))**2)) * np.sqrt(np.sum((y-np.mean(y))**2)))

    if not CI:
        return rho

    z = norm.ppf(sum([CI, 1.0])/2)
    lower = np.tanh(np.arctanh(rho) - z/np.sqrt(n-3))
    upper = np.tanh(np.arctanh(rho) + z/np.sqrt(n-3))
    return (rho, lower, upper)
