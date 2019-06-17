import numpy as np

from semiconductor_photophysics import ini_parsing

# --- helper methods to convert from unit-ed system to unitless system ----------------------------

kb = 8.6173e-5  # ev/K
kb_J = 1.381e-23  # J/K
hbar_Js = 1.055e-34  # Js
m_e = 9.1094e-31  # kg


def n_cm3_to_nm3(n_cm):
    return n_cm / (1e7) ** 3  # nm^-3


def E_to_wbar(E, Eg, ER):
    return (E - Eg) / ER


def T_to_Tbar(T, ER):
    return T * kb / ER


def G_to_Gbar(G, ER):
    return G / ER


def dielectric_prefactor(rcv, a0, ER):
    return rcv ** 2 / (2 * np.pi * a0 ** 3 * ER)


def mubar(mu, Eg, T):
    return (mu - Eg / 2) / (kb * T)


def g_from_ak(a0, k):
    return 12 / (np.pi ** 2 * a0 * k)


def g_from_exp_params(n, epsilonr, T, m_star):
    a = a0_exciton(epsilonr, m_star)
    k = k_debye_huckel(n, epsilonr, T)
    return g_from_ak(a, k)


def Eg_from_g(Eg0, ER, g, bound=True):
    if bound:
        out = Eg0 + ER * (-1 + (1 - 1 / g) ** 2)
    else:
        out = Eg0 - ER / g
    return out


def k_debye_huckel(n, epsilonr, T):
    # n has units of 1/nm^3
    hc = 1240.  # ev nm
    beta = 1 / (kb * T)  # 1/eV
    alpha = 1 / 137.036
    return np.sqrt(2 * n * beta * alpha * hc / epsilonr)


def a0_exciton(epsilonr, m_star):
    # epsilonr is unitless
    # m_star is unitless reduced mass of exciton hole system m_star = m/me
    a0 = 5.29e-2  # nm
    return a0 * epsilonr / m_star


def ER_exciton(a0_exc, m_star):
    # a0_exc is excition bohr radius in nm
    # m_star is unitless reduced mass of exciton hole system m_star = m/me
    m = m_star * m_e  # kg
    a0 = a0_exc * 1e-9  # nm to m
    out = hbar_Js ** 2 / (2 * m * a0 ** 2)  # J
    return out * 6.242e18  # eV


def m_star_calc(me_star, mh_star):
    return me_star * mh_star / (me_star + mh_star)


def n0_nm(T, m_alpha_star):
    m = m_alpha_star * m_e  # kg
    beta = 1 / (kb_J * T)
    out = .25 * (2 * m / (hbar_Js ** 2 * np.pi * beta)) ** (3 / 2)  # m^-3
    out /= (1e9) ** 3  # nm^-3
    return out


def mu_from_n(n, m_alpha_star, T):
    n0 = n0_nm(T, m_alpha_star)
    nu = n / n0
    # from Haug and Koch pg 96
    K1 = 4.897
    K2 = 0.045
    K3 = 0.133
    mu = kb * T * (np.log(nu) + K1 * np.log(K2 * nu + 1) + K3 * nu)
    return mu  # eV


# --- methods for calculating dielectric spectrum according to Banyai and Koch --------------------


def L(E, E0, G, A=1., A0=0.):
    """
    Complex Lorentzian. Area normalized to imaginary component.
    """
    return A / np.pi / (E0 - E - 1j * G) + A0


def _checkndim_copy_reshape(arrs, added_dims):
    arrs_out = []
    for arr in arrs:
        assert (
            arr.ndim == arrs[0].ndim
        ), "Inputs need to be arrays all with the same ndim"
        arr_out = arr.copy()
        arr_out = np.reshape(arr_out, arr.shape + (1,) * added_dims)
        arrs_out.append(arr_out)
    return arrs_out


def band_filling_factor(w, T, mu_e, mu_h):
    argument = w / T - mu_e - mu_h
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
    #* L(w + (1 / ell - ell / g) ** 2, 0, G)
    def f2(w, g, G, ell):
        return (
            np.pi
            * L(w + (1 / ell - ell / g) ** 2, 0,G)
            * 2
            * (g - ell ** 2)
            * (2 * ell ** 2 - g)
            / (ell ** 3 * g ** 2)
        )

    def f3(ell):
        out = np.ones(ell.shape)
        out[ell == 0] = 0
        return out

    def f4(g, ell, n):
        # we aren't ensuring that we don't compute ell == n.
        # TODO catch ell == n with a more robust method
        out = (
            n ** 2
            * (n ** 2 * ell ** 2 - (g - ell ** 2) ** 2)
            / ((n ** 2 - ell ** 2) * (n ** 2 * ell ** 2 - g ** 2))
        )
        out[out == -np.inf] = np.nan
        out[out == np.inf] = np.nan
        return out

    # check input arrays
    ndim_org = w.ndim
    w, g, G = _checkndim_copy_reshape([w, g, G], 2)
    # create array to sum over
    ellmax = int(np.max(np.sqrt(g)))
    ell_base = np.arange(1, ellmax + 1, 1)
    new_shape = (1,) * ndim_org + (ell_base.size,) + (1,)
    ell_base = np.reshape(ell_base, new_shape)
    ell = np.repeat(
        ell_base, g.size, axis=1
    )  # TODO np.argmax g.shape | assumes g is varying along 1st axis
    ell[f1(g) < ell] = 0
    # create array to product over
    n = np.arange(1, nmax + 1, 1)
    new_shape = (1,) * (ndim_org + 1) + n.shape
    n = np.reshape(n, new_shape)
    # now actually create output arrays
    out = f4(g, ell, n)
    out = np.nanprod(out, axis=-1, keepdims=True)
    out = out * f2(w, g, G, ell) * f3(ell)
    out = np.nansum(
        out, axis=-2, keepdims=True
    )  # TODO figure out why I need nansum here
    if squeeze:
        out = np.squeeze(out)
    return out


def continuum_contribution(w, g, G, xmax, xnum, nmax, squeeze=True):
    # TODO docstring
    """
    """
    w, g, G = _checkndim_copy_reshape([w, g, G], 2)
    # create n array to take product against and x array to integrate against
    # arrays have shape (..., x, n) where ... is specified by w,g,G arrays.
    x = np.linspace(0, xmax, xnum)
    dx = x[1] - x[0]
    x = np.reshape(x, (1,) * w.ndim + x.shape + (1,))
    n = np.arange(1, nmax + 1, 1)
    n = np.reshape(n, (1,) * (w.ndim + 1) + n.shape)
    # Coulomb enhancement product
    coul_enhnc = 1 + (2 * g * n ** 2 - g ** 2) / (
        (n ** 2 - g) ** 2 + n ** 2 * g ** 2 * x
    )
    coul_enhnc = np.prod(coul_enhnc, axis=-1, keepdims=True)
    # create array to integrate
    arr = np.sqrt(x) * coul_enhnc * L(w-x, 0, G) # different than BK, but they are working with delta functions we care about the sign
    out = np.trapz(arr, dx=dx, axis=-2)
    if squeeze:
        out = np.squeeze(out)
    return out


def reduced_dielectric(wbar, g, Gbar, Tbar, mubar_e, mubar_h, xmax, xnum, nmax):
    # TODO docstring
    """
    reduced units. no prefactor
    """
    bf = band_filling_factor(wbar, Tbar, mubar_e, mubar_h)
    bound = bound_contribution(wbar, g, Gbar, nmax, squeeze=True)
    continuum = continuum_contribution(wbar, g, Gbar, xmax, xnum, nmax, squeeze=True)
    out =  bf * (bound + continuum)
    #print(bound.imag.max(), continuum.imag.max())
    return out

def dielectric_microscopic(w, Eg0, G, a0, k, T, rcv, mu_e, mu_h, m_star, xmax, xnum, nmax, print_g=False, return_bff=False):    
    # TODO docstring
    """
    """
    g = g_from_ak(a0, k)
    if print_g:
        print("g min and max", g.min(), g.max())
    ER = ER_exciton(a0, m_star)
    if g.min() < 1:
        raise Exception("currently g less than 1 is not implemented---this is a Mott transition")
    Eg = Eg_from_g(Eg0, ER, g, bound=True)
    Gbar = G_to_Gbar(G, ER)
    Tbar = T_to_Tbar(T, ER)
    wbar = E_to_wbar(w, Eg, ER)
    pre = dielectric_prefactor(rcv, a0, ER)
    mubar_e = mubar(mu_e, Eg, T)
    mubar_h = mubar(mu_h, Eg, T)
    out = pre * reduced_dielectric(wbar, g, Gbar, Tbar, mubar_e, mubar_h, xmax, xnum, nmax)
    bf = band_filling_factor(wbar, Tbar, mubar_e, mubar_h)
    if return_bff:
        return out, bf
    else:
        return out

def dielectric_microscopic_from_ini(w, p, out_shape=None):
    params = ini_parsing.read_full_sim_params(p)
    num_params = params['num_params']
    if out_shape == None:
        out = np.zeros(w.shape, dtype=complex)
    else:
        out = np.zeros(out_shape, dtype=complex)
    for BK in params['BKs']:
        out += dielectric_microscopic(w, *BK, *num_params, return_bff=False)
    for lor in params['Lors']:
        out += L(w, *lor)
    return out  

def dielectric_macroscopic():
    # TODO docstring
    """
    """
    pass

