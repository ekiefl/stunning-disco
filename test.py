from colour import Color as colour

def get_color_gradient(base_colors, gradations=100):
    """
    This function returns a list of all the colors for a multi-color gradient.
    
    INPUTS
    ------
    base_colors : list
        This is a list of strings interpretable by the `colour` module as colors.  
        If you want a 3-color gradient from red to white to blue, then base_colors
        should be ["red", "white", "blue"], or [#ff0000, #ffffff, #0000ff].
    gradations : integer
        The approximate number of gradations asked for (It will not be exactly this
        unless the modulus of gradiations%(n-1) == 0.)
    """
    
    n = len(base_colors)

    if gradations // (n-1) < 2:
    #   there is no work to do. just return the colors in base_colors as Color objects
        return [colour(base_color) for base_color in base_colors]

    color_gradient = []

    for i in range(1, n):

        fro = colour(base_colors[i-1])
        to = colour(base_colors[i])

        color_gradient_interval = list(fro.range_to(to, gradations//(n-1)))   
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

base_colors = ["red","blue","green"]
color_gradient = get_color_gradient(base_colors)
print(len(color_gradient))




