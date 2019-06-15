""" Elliott model of absorption. """
# TODO add references

# --- import ------------------------------------------------------------------

import numpy as np

import scipy.signal as sig


# --- define ------------------------------------------------------------------


__all__ = ['meh']


# --- Absorption class ----------------------------------------------------


class Absorption():
    
    __init__(self, Material):
        

    
    def c_e(self, D):
        """
        Calculate Coulomb enhancement factor. 
        """
        D[D<0] = 0.0 # coerce all negative inputs to zero 
        d = np.pi/np.sqrt(D)
        D[np.isnan(D)] = 0.
        return np.pi*np.exp(d)/np.sinh(d)



    def exciton_series(self, D, FWHM, j):
        """
        Generates Gaussian resonances from hydrogenic seres of cardinality j
        
        D is in normalized and shifted energy units. FWHM is not. 
        """
        D = check_arr(D)
        FWHM = check_arr(FWHM)
        js = np.linspace(1,j,j)
        out = 4*np.pi/js[None,None,:]**3 * gauss(D[:,None,None]+1/js[None,None,:]**2,FWHM[None,:,None])
        out = np.sum(out, axis=-1)
        return out # 2D array

    def continuum(self, D, FWHM_cont):
        """
        Generates continuum contribution to spectra
        """
        d = np.pi/np.sqrt(D)
        out = np.heaviside(D,.5) * np.pi * np.exp(d) / np.sinh(d)
        out[np.isnan(out)] = 0.
        #original_max = out.max()
        g = gauss(D-D.mean(), FWHM_cont)
        g /= g.sum() # This normalization is needed so that the convolution is invarient to number of x points
        out = sig.fftconvolve(g, out, mode='same')
        #out /= original_max
        return out

    def alpha(self, E, Eg, Ebind, FWHM_ex, j, FWHM_cont=0, A=1, B=1):
        D = (E-Eg)/Ebind
        ex = exciton_series(D, FWHM_ex/Ebind, j)
        ex = ex[...,-1]
        contium = continuum(D, FWHM_cont/Ebind)
        return E/Ebind * (A*ex + B*contium)