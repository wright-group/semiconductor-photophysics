""" Optical transforms. """

# --- import ------------------------------------------------------------------

import numpy as np

import scipy.signal as sig


# --- define ------------------------------------------------------------------


__all__ = ['KK']


# --- Kramers-Kronig class ----------------------------------------------------


class KK():
    
    def __init__(self, alpha, E, ell=1e-5, axis=0, offset=0):
        self.alpha = alpha # in units of inverse cm
        self.E = E # in units of eV
        self.ell = ell # sample length in units of cm
        self.hbarc = 1.9732697e-05 # eV cm
        self.axis = axis
        self.k_i()
        self.offset = offset
        self.n_r()
        self.R_normal()
        self.T_normal()
        
    def k_i(self):
        """
        Calculate absorptive part of refractive index.
        """
        self.k = self.hbarc/2/self.E * self.alpha

    def n_r(self):
        """
        Calculate dispersive part of refractive index using Hilbert transform.
        """
        dE = self.E[1]-self.E[0] # step size
        # generate mirrored and flipped imaged array because k is globaly odd
        k_full = np.concatenate((-np.flip(self.k,axis=self.axis), self.k), axis=self.axis) 
        # normalize so that the Hilbert transform is invarient to number of x points
        k_full /= self.k.sum(axis=self.axis)[None,...]
        # Do Hilbert transform
        n = 1 - np.imag(sig.hilbert(k_full, axis=self.axis))/dE + self.offset
        self.n = n[self.E.size :,...] # dump the negative frequency image

    def R_normal(self):
        """
        Calculate reflection with beam normal to surface.
        """
        self.R = ((self.n-1)**2 + self.k**2) / ((self.n+1)**2 + self.k**2)
    
    def T_normal(self):
        """
        Calculate transmission with beam normal to surface.
        
        See Pankove's 'Optical Processes in Semiconductors' eq 4-32.
        """
        self.T = (1-self.R)**2 * np.exp(-self.alpha * self.ell) / (1-self.R**2*np.exp(-2*self.alpha*self.ell))