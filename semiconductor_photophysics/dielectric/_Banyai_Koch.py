import numpy as np
import matplotlib.pyplot as plt

def L(E, E0, G):
    """
    Complex Lorentzian. Area normalized to imaginary component.
    """
    return 1/np.pi / (E0-E-1j*G)


def _checkndim_copy_reshape(arrs, added_dims):
    arrs_out = []
    for arr in arrs:
        assert arr.ndim == arrs[0].ndim, "Inputs need to be arrays all with the same ndim"
        arr_out = arr.copy()
        arr_out = np.reshape(arr_out, arr.shape + (1,)*added_dims)
        arrs_out.append(arr_out)
    return arrs_out



def band_falling_factor(w, T, mu_e, mu_h):
    argument = w/T - mu_e - mu_h
    argument *= .5
    return np.tanh(argument)


def bound_contribution(w, g, G, nmax, squeeze=True):
    # TODO docstring
    # This method assumes g, if it is multidimensional varies along the first axis while w varies along the zeroth axis.
    """
    """
    # helper methods
    def f1(g):
        return np.sqrt(g).astype(int)    
    def f2(w, g, G, ell):
        return np.pi * L( w+(1/ell - 1/g)**2, 0, G) * 2*(g-ell**2)*(2*ell**2-g)/(ell**3*g**2)    
    def f3(ell):
        out = np.ones(ell.shape)
        out[ell == 0] = 0
        return out    
    def f4(g, ell, n):
        # we aren't ensuring that we don't compute ell == n. 
        # TODO catch ell == n with a more robust method 
        out = n**2 * (n**2*ell**2 - (g-ell**2)**2) / ((n**2-ell**2) * (n**2*ell**2 - g**2))
        out[out == -np.inf] = np.nan
        out[out == np.inf] = np.nan
        return out
    # check input arrays
    ndim_org = w.ndim
    w, g, G = _checkndim_copy_reshape([w,g,G], 2)
    # create array to sum over
    ellmax = int(np.max(np.sqrt(g)))
    ell_base = np.arange(1,ellmax+1,1)
    new_shape = (1,)*ndim_org + (ell_base.size,) + (1,)
    ell_base = np.reshape(ell_base, new_shape)
    ell = np.repeat(ell_base, g.size, axis=1) # TODO np.argmax g.shape | assumes g is varying along 1st axis
    ell[f1(g) < ell] = 0
    # create array to product over
    n = np.arange(1, nmax + 1, 1)
    new_shape = (1,)*(ndim_org + 1) + n.shape
    n = np.reshape(n, new_shape)   
    # now actually create output arrays
    out = f4(g, ell, n)
    out = np.nanprod(out, axis=-1, keepdims=True)
    out = out * f2(w, g, G, ell) * f3(ell)
    out = np.nansum(out, axis=-2, keepdims=True) #TODO figure out why I need nansum here
    if squeeze:
        out = np.squeeze(out)
    return out       


def continuum_contribution(w, g, G, xmax, xnum, nmax, squeeze=True):
    # TODO docstring
    """
    """
    w, g, G = _checkndim_copy_reshape([w,g,G], 2) 
    # create n array to take product against and x array to integrate against
    # arrays have shape (..., x, n) where ... is specified by w,g,G arrays.
    x = np.linspace(0,xmax,xnum)
    dx = x[1] - x[0]
    x = np.reshape(x, (1,)*w.ndim + x.shape + (1,))
    n = np.arange(1, nmax + 1, 1)
    n = np.reshape(n, (1,)*(w.ndim + 1) + n.shape)  
    # Coulomb enhancement product
    coul_enhnc = 1 + (2*g*n**2 -g**2) / ((n**2 - g)**2 + n**2*g**2*x)
    coul_enhnc = np.prod(coul_enhnc, axis=-1, keepdims=True)
    # create array to integrate
    arr = np.sqrt(x) * coul_enhnc * L(x-w, 0, G)
    out = np.trapz(arr, dx=dx, axis=-2)
    if squeeze:
        out = np.squeeze(out)
    return out
   
    



    
    