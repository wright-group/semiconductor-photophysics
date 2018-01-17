"""Tanguy model of dielectric function. """
# DOI: 10.1103/PhysRevLett.75.4090

# --- import --------------------------------------------------------------------------------------

import numpy as np

import matplotlib.pyplot as plt

from scipy.special import digamma as psi


# --- define --------------------------------------------------------------------------------------


__all__ = ['meh']


# --- workspace -----------------------------------------------------------------------------------


def xi(z, Eg, R):
    """
    """
    return np.sqrt(R/(Eg-z))

def g(z, Eg, R, allowed=True):
    """
    """
    xi_ = xi(z, Eg, R) # only evalute once
    ga = 2*np.log(xi_) - 2*np.pi*(np.cos(np.pi*xi_)/np.sin(np.pi*xi_)) - 2*psi(xi_) - 1/xi_
    if allowed:
        return ga
    else:
        gf = (z - Eg + R)*ga
        return gf

def e_G(E, Eg, R, G, A, allowed=True):
    """
    """
    z = E + 1j*G
    prefactor = A * np.sqrt(R) / (z)**2
    postfactor = g(z, Eg, R, allowed=allowed) 
    postfactor += g(-1*z, Eg, R, allowed=allowed) 
    postfactor -= 2 * g(0, Eg, R, allowed=allowed)
    return prefactor * postfactor

def e_both(E, Eg, R, G, A_allowed, A_forbidden):
    """
    """
    return e_G(E, Eg, R, G, A_allowed, allowed=True) + e_G(E, Eg, R, G, A_forbidden, allowed=False)

def complex_index(e):
    """
    """
    mag = np.abs(e)
    e1 = np.real(e)
    n = np.sqrt((mag + e1)/2)
    k = np.sqrt((mag - e1)/2)
    return n, k

x = np.linspace(1.5,2.5,1000)
y = e_both(x, 2.4, .5, .05, 4, 2) + e_both(x, 2.4, .4, .05, 4, 2)
n, k = complex_index(y)
#y = e_G(x, 2, .1, .01, 1, allowed=False)
#plt.plot(x, y.real)
#plt.plot(x, y.imag)

plt.plot(x, n)
plt.plot(x, k)