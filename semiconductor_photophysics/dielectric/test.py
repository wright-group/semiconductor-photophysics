import _Banyai_Koch as BK
import numpy as np
import WrightTools as wt
import matplotlib.pyplot as plt

# --- literature values for MAPbI3 perovskites ----------------------------------------------------

m_star = 0.104 # 10.1038/NPHYS3357
ER = .02#.003 # 10.1038/NPHYS3357 reports at value at RT as "a few meV"
epsilonr = 10#25.7 # 10.1103/PhysRevB.89.155204 theory paper. 0 frequency dielectric. Other papers say mid frequency dielectric should be used instead
Eg0 = 1.6
a0 = BK.a0_exciton(epsilonr, m_star)

#ER = BK.ER_exciton(a0, m_star)

# --- throw-away values




w = np.linspace(1.,3,1001)
G = np.array([.01])
T = 300
rcv = 1.
mu_e = .7
mu_h = .7
xmax = 1000
xnum = 10000
nmax = 1000

n_nm = BK.n_cm3_to_nm3(1e17)
k = BK.k_debye_huckel(n_nm, epsilonr, T)
k = 1/100.
g = BK.g_from_ak(a0, k)

a0 = np.array([a0])
k = np.array([k])

if True:
    out = BK.dielectric_microscopic(w, Eg0, G, a0, k, T, rcv, mu_e, mu_h, m_star, xmax, xnum, nmax)
    plt.plot(w, out.imag)
    plt.plot(w, out.real)
    plt.plot(w, np.abs(out))

if False:
    Tbar = BK.T_to_Tbar(300, ER)
    wbar = BK.E_to_wbar(w, Eg0, ER)
    mubar_e = BK.mubar(mu_e, Eg0, T)
    mubar_h = BK.mubar(mu_h, Eg0, T)    
    out = BK.band_filling_factor(wbar, Tbar, mubar_e, mubar_h)
    plt.plot(w,out)
    
