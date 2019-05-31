""" 
Implementation of Matrix Transfer approach to calculating transmission, reflection, and absorption.
Many ideas were taken from Steven Byrnes implementation in the tmm package 
https://github.com/sbyrnes321/tmm/blob/master/tmm_core.py (used under terms of MIT license)
The current implementation allows for multidimensional dielectric arrays.

"""

import numpy as np
import scipy as sp



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

def _r_from_M():
    return M[...,1,0] / M[...,0,0]

def _t_from_M():
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

def _M_generator(pol, n_i, n_f, th_i, th_f, delta):
    # eq 11 in byrnes notes
    rnn1 =  _r_calc(pol, n_i, n_f, th_i, th_f)
    tnn1 =  _t_calc(pol, n_i, n_f, th_i, th_f)
        
def _snells_law_calc(n_1, n_2, th_1):
    # TODO this is super naive. Consider making sure we are in the correct branch cut
    th_2_guess = sp.arcsin(n_1*np.sin(th_1) / n_2)
    return th_2_guess

def stack_calculation(pol, n_arr, d_arr, th_0, hw_vac):
    # TODO consider switching to lists of arrays
    lam_vac = wt.units.converter(hw_vac, 'eV', 'nm')
    
    # TODO calculate angles
    th_arr 
    
    kz_arr = 2 * np.pi * n_arr * np.cos(th_arr) / lam_vac
    delta_arr = kz_arr * d_arr
    
    
    Mlist = []
    M_0 =
    Mlist.append(M0)
    # TODO build rest of M
    
    
    
    Mout = _Mlist_prod(Mlist)
    