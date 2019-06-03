""" 
Implementation of Matrix Transfer approach to calculating transmission, reflection, and absorption.
Many ideas were taken from Steven Byrnes implementation in the tmm package 
https://github.com/sbyrnes321/tmm/blob/master/tmm_core.py (used under terms of MIT license)
The current implementation allows for multidimensional dielectric arrays.

"""

import numpy as np
import scipy as sp
import WrightTools as wt



def _R_from_r(r):
    """
    Calculate reflected power R, starting with reflection amplitude r.
    """
    return np.abs(r)**2

def _T_from_t(pol, t, n_i, n_f, th_i, th_f):
    """
    Calculate transmitted power T, starting with transmission amplitude t.
    """
    if pol == 's':
        return np.abs(t**2) * (((n_f*np.cos(th_f)).real) / (n_i*np.cos(th_i)).real)
    elif pol == 'p':
        return np.abs(t**2) * (((n_f*np.conj(np.cos(th_f))).real) / (n_i*np.conj(np.cos(th_i))).real)
    else:
        raise ValueError("Polarization must be 's' or 'p'")    

def _r_from_M(M):
    return M[...,1,0] / M[...,0,0]

def _t_from_M(M):
    return 1 / M[...,0,0]

def _Mlist_prod(Mlist):
    Mout = Mlist.pop(0)
    for M in Mlist:
        Mout = Mout @ M
    return Mout

def _t_calc(pol, n_i, n_f, th_i, th_f):
    if pol == 's':
        return 2 * n_i * np.cos(th_i)  / ( n_i * np.cos(th_i) + n_f * np.cos(th_f) )
    elif pol == 'p':
        return 2 * n_i * np.cos(th_i)  / ( n_f * np.cos(th_i) + n_i * np.cos(th_f) )
    else:
        raise ValueError("Polarization must be 's' or 'p'")  

def _r_calc(pol, n_i, n_f, th_i, th_f):
    if pol == 's':
        out = n_i * np.cos(th_i) - n_f * np.cos(th_f)
        out /= n_i * np.cos(th_i) + n_f * np.cos(th_f)
        return out
    elif pol == 'p':
        out = n_f * np.cos(th_i) - n_i * np.cos(th_f)
        out /= n_f * np.cos(th_i) + n_i * np.cos(th_f)
        return out
    else:
        raise ValueError("Polarization must be 's' or 'p'")       

def _M_generator(pol, n_i, n_f, th_i, th_f, deltan):
    # eq 11 in byrnes notes
    rnn1 =  _r_calc(pol, n_i, n_f, th_i, th_f)
    tnn1 =  _t_calc(pol, n_i, n_f, th_i, th_f)

    M1 = np.zeros(deltan.shape + (2,2), dtype=complex)
    M1[...,0,0] = np.exp(-1j*deltan) # TODO ensure matrix construction is as intended 
    M1[...,1,1] = np.exp(1j*deltan) # TODO ensure matrix construction is as intended 
    
    M2 = np.ones(deltan.shape + (2,2), dtype=complex)
    M2[...,0,1] = rnn1 # TODO ensure matrix construction is as intended 
    M2[...,1,0] = rnn1 # TODO ensure matrix construction is as intended 

    out = M1 @ M2
    out /= tnn1[..., None, None]
    return out

def _M_bootstrap(pol, n, th, deltan):
    assert n.shape == th.shape == deltan.shape, "input arrays have mismatched shapes" 
    Mout = []
    for i in range(1, n.shape[0]):
        M = _M_generator(pol, n[i-1], n[i], th[i-1], th[i], deltan[i-1])
        Mout.append(M)
    return Mout
    
        
def _snells_law_calc(n_1, n_2, th_1):
    # TODO this is super naive. Consider making sure we are in the correct branch cut
    th_2_guess = sp.arcsin(n_1*np.sin(th_1) / n_2)
    return th_2_guess

def _snells_bootstrap(ns, th_0):
    theta_out = np.zeros(ns.shape, dtype=complex)
    theta_out[0] = th_0
    for i in range(1, ns.shape[0]):
        theta_old = theta_out[i-1]
        n_old = ns[i-1]
        n_new = ns[i]
        theta_new = _snells_law_calc(n_old, n_new, theta_old)
        theta_out[i] = theta_new
    return theta_out

def stack_calculation(pol, n_arr, d_arr, th_0, hw_vac):
    """ Calculate optical properties of a stack of optical structures.
    This calculator assumes arrays are well shaped. 
    0th dimension of arrays correlate to optical stack number.
    1st dimension of arrays correlate to energy/wavelength of light
    2nd and more dimensions correlate to user specified refractive index changes
    
    Parameters
    ----------
    pol : string
        's' or 'p' specifies the polarization type
    n_arr : array
        refractive indecies of optical stack
        For x layers (include the input, leading, and output, trailing, layers) required shape is
        (x, y, ...).
        By convention, the first and last layers have exclusively real refractive indecies.
    d_arr : array
        thicknesses of optical stack in nanometers.
        For x layers required shape is (x, 1, ...).
        By convention, first and last layers have zero thickness. 
    th_0 : float
        angle of forward traveling light from 0th to 1st layer
    hw_vac : array
        energy per photon of light in vacuum
        must be of shape (1, y, ...)
        
    Returns
    -------
    tuple
        R, T, A: arrays
            arrays have shape (y, ...)
            R : reflectance
            T : transmittance
            A : absorptance
    
    """
    # ensure d_arr has zero thickness for first and last layers
    d_arr[0] = 0
    d_arr[-1] = 0
    # convert to nm
    lam_vac = wt.units.converter(hw_vac, 'eV', 'nm')
    # calculate arrays
    th_arr = _snells_bootstrap(n_arr, th_0)    
    kz_arr = 2 * np.pi * n_arr * np.cos(th_arr) / lam_vac
    delta_arr = kz_arr * d_arr
    # create list of M arrays
    Mlist = _M_bootstrap(pol, n_arr, th_arr, delta_arr)
    # now take their product    
    Mout = _Mlist_prod(Mlist)
    # calculate useful quantities
    r = _r_from_M(Mout)
    t = _t_from_M(Mout)
    R = _R_from_r(r)
    T = _T_from_t(pol, t, n_arr[0], n_arr[-1], th_arr[0], th_arr[-1])
    A = 1 - R - T
    return R, T, A
 
    

# TODO remove
if False:
    import matplotlib.pyplot as plt
    w = np.linspace(1,3, 301)[:,None] * np.ones((301,201))
    w1 = np.linspace(1,3,201)[None,:] * np.ones((301,201))
    zero = np.ones(w.shape)
    first = 1.5 + .01 / (w-2-1j*.05)
    first2 = 1.5 + .01 / (w-2-1j*.05) * w1
    second = np.ones(w.shape)
    
    # x, 31, 21
    
    arrs = (zero, first, first2, second)
    narr = np.stack(arrs)
    d_arr = np.array([0, 100, 100, 0])[:,None, None]
    R, T, A = stack_calculation('s', narr, d_arr, 0.5, w)
    
    #plt.plot(w, R)
    #plt.plot(w, T)
    #plt.plot(w, A)
    
    
    