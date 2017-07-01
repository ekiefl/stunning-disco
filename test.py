from colour import Color as colour
import numpy as np
import matplotlib.pyplot as plt
import pylab

def get_color_gradient(base_colors, gradations=100):
    """
    This function returns a list of all the colors for a multi-color gradient.
    
    INPUTS
    ------
    base_colors : list
        This is a list of strings interpretable by the `colour` module as colors.  
        If you want a 3-color gradient from red to white to blue, then base_colors
        should be ["red", "white", "blue"], or [#ff0000, #ffffff, #0000ff].
    gradations : integer, list
        The approximate number of gradations asked for. If an integer is provided, 
        a color list with len(gradations) will be returned (It will not be exactly this
        unless the modulus of gradiations%(n-1) == 0). If you want the number of
        gradations to vary between the colors in base_colors, you can provide gradations
        as a list. In this case, the nth element in gradations corresponds to the number
        of gradations in the color range defined by the nth and (n+1)th colors in 
        base_colors. For example, if base_colors = ["red","blue","green"] and base_colors
        = [10, 5], 10 colors are defined between red and blue and 5 colors between blue and
        green. The sum of the entries are therefore the total number of color gradations
        expected.
    """

    
    n = len(base_colors)

#   if int-like gradations is passed, transform into array-like
    if type(gradations) == int:
        gradations = [gradations//(n-1) for _ in range(n-1)]

#   if array-like gradations is passed, ensure its length is n-1
    if len(gradations) != n-1:
            raise ValueError("The gradations array must be n-1")

    if any([x==0 for x in gradations]):
        raise ValueError("One or more of the gradation intervals you defined sucks.")

    color_gradient = []

    for i in range(1, n):

        fro = colour(base_colors[i-1])
        to = colour(base_colors[i])

        color_gradient_interval = list(fro.range_to(to, gradations[i-1]))
        if not i == n - 1: del color_gradient_interval[-1]
        color_gradient.extend(color_gradient_interval)

    return color_gradient

def map(x, color_gradient, m=0.0, M=1.0):
    """ 
    This map converts a numeric value to Color object in a list of Color objects 
    """
    n = len(color_gradient)
    numeric_range = np.linspace(m,M,n)
    nearest_idx = np.abs(x - numeric_range).argmin()
    return color_gradient[nearest_idx]

base_colors = ["#a03129","#fcf5f4","#ffffff","#e8edf9", "#264799"]
gradations = [50,25,25,50]

color_gradient = get_color_gradient(base_colors, gradations)

print(color_gradient)
print(len(color_gradient))

x = np.arange(len(color_gradient))
pylab.show()
plt.scatter(x,x,color=[x.rgb for x in color_gradient])



