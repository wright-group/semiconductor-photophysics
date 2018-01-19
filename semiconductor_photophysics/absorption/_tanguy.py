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

def g(z, Eg, R, allowed=True, three_D=True):
    """
    """
    xi_ = xi(z, Eg, R) # only evalute once
    if three_D:
        ga = 2*np.log(xi_) - 2*np.pi*(np.cos(np.pi*xi_)/np.sin(np.pi*xi_)) - 2*psi(xi_) - 1/xi_
        if allowed:
            return ga
        else:
            gf = (z - Eg + R)*ga
            return gf
    else:
        ga = 2*np.log(xi_) - 2*psi(0.5-xi_)
        if allowed:
            return ga
        else:
            gf = (z - Eg + 4*R)*ga
            return gf
        
def e_G(E, Eg, R, G, A, allowed=True, three_D=True):
    """
    """
    z = E + 1j*G
    if three_D:
        prefactor = A * np.sqrt(R) / (z)**2
    else:
        prefactor = A  / (z)**2 / np.pi
    postfactor = g(z, Eg, R, allowed=allowed, three_D=three_D) 
    postfactor += g(-1*z, Eg, R, allowed=allowed, three_D=three_D) 
    postfactor -= 2 * g(0, Eg, R, allowed=allowed, three_D=three_D)
    return prefactor * postfactor

def e_both(E, Eg, R, G, A_allowed, A_forbidden, three_D=True):
    """
    """
    return e_G(E, Eg, R, G, A_allowed, allowed=True, three_D=three_D) + e_G(E, Eg, R, G, A_forbidden, allowed=False, three_D=three_D)

def complex_index(e):
    """
    """
    mag = np.abs(e)
    e1 = np.real(e)
    n = np.sqrt((mag + e1)/2)
    k = np.sqrt((mag - e1)/2)
    return n, k

if False:
    x = np.linspace(1.5,2.5,1000)
    y3 = e_both(x, 2.4, .5, .05, 4, 2, three_D=True) + e_both(x, 2.4, .4, .05, 4, 2, three_D=True)
    y2 = e_both(x, 2.4, .5, .05, 4, 2, three_D=False)# + e_both(x, 2.4, .4, .05, 4, 2, three_D=False)
    # need to test 2D. Reproduce results of  Tanguy's paper.
    plt.plot(x, y3.real)
    plt.plot(x, y3.imag)
    
    #plt.plot(x, y2.real)
    plt.plot(x, y2.imag)
    
if False:
    x = np.linspace(1.35,1.45, 1000)
    y = e_both(x, 1.42, 0.004, 0.006, 1, 0, False)
    #y /= np.sqrt(0.004)
    reduced = (x-1.42)/.004
    plt.plot(reduced, y.real)
    plt.plot(reduced, y.imag)
    
x = np.linspace(1,2, 1000)
y1 = -1*np.exp(1j*0)*(x-1.5+1j*0.05)**(-.5)
y2 = 3-1*np.exp(1j*0)*np.log(x-1.5+1j*0.05)
y3 = -1*np.exp(1j*0)*(x-1.5+1j*0.05)**(.5)

plt.plot(x, y1.imag, color='C0')
plt.plot(x, y1.real, color='C0')
plt.plot(x, y2.imag, color='C1')
plt.plot(x, y2.real, color='C1')
plt.plot(x, y3.imag, color='C2')
plt.plot(x, y3.real, color='C2')