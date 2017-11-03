""" Toolbox of useful things """

# --- import ------------------------------------------------------------------


# --- constants ---------------------------------------------------------------


k = 8.6173303e-5 # Boltzmann constant in units of eV/K


# --- helper methods ----------------------------------------------------------


def check_arr(i):
    """
    Coerces input to array type.
    """
    if not hasattr(i, 'size'):
        i = np.array([i], dtype=float)
    return i


def clipper(arr, frac=.01):
    """
    Clips first dimension of array by some fraction of the length of the dimension
    """
    tot = arr.shape[0]    
    num = int(tot*frac) 
    return arr[:-num,...]


def gauss(x, FWHM, x0=0):
    """
    Returns Gaussian lineshape. 
    """
    s = FWHM/2.355
    return np.exp(-(x-x0)**2/(2*s**2))


def lorentz(x, FWHM, x0=0):
    """
    Returns Lorentzian lineshape.
    """
    return FWHM/2/np.pi*((x-x0)**2+(FWHM/2)**2)**-1


    
