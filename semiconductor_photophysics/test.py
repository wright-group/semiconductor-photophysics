import dielectric_Banyai_Koch as BK
import optics as optics
import numpy as np
import WrightTools as wt
import matplotlib.pyplot as plt

# --- literature values for MAPbI3 perovskites ----------------------------------------------------

m_star = 0.104 # 10.1038/NPHYS3357
ER = .02#.003 # 10.1038/NPHYS3357 reports at value at RT as "a few meV"
epsilonr = 9#25.7 # 10.1103/PhysRevB.89.155204 theory paper. 0 frequency dielectric. Other papers say mid frequency dielectric should be used instead
Eg0 = 1.67
a0 = BK.a0_exciton(epsilonr, m_star)
#a0 = 5

#ER = BK.ER_exciton(a0, m_star)

# --- throw-away values




w = np.linspace(1.3,3.5,1501)
G = np.array([.04])
T = 300
rcv = 1.75
mu_e = .1
mu_h = .1
xmax = 200
xnum = 1000
nmax = 100

n_nm = BK.n_cm3_to_nm3(1e15)
k = BK.k_debye_huckel(n_nm, epsilonr, T)
#k = 1/100.
g = BK.g_from_ak(a0, k)

a0 = np.array([a0])
k = np.array([k])

if True:
    #col = wt.open('perovskite_abs.wt5')
    #d = col.abs_film1
    #values = d.signal[:]
    #values -= values.min()
    #values = 1 - 10**(-1*values)
    #d.create_channel(name='frac_abs', values=values)
    
    
    
    wreshape = w[None,:]
    out = 1.8*BK.dielectric_microscopic(w, Eg0, G, a0, k, T, rcv, mu_e, mu_h, m_star, xmax, xnum, nmax) 
    out += BK.L(w, 2.15, .2)*.65
    out += BK.L(w, 2.5, .2)*.5
    out += BK.L(w, 3.3, .2)*1.5
    out += 10
    n_out = optics.e_to_n(out)
    #w = np.linspace(1,3, 301)[:,None] * np.ones((301,201))
    #w1 = np.linspace(1,3,201)[None,:] * np.ones((301,201))
    zero = 1.52*np.ones(w.shape)
    first = n_out
    second = 1.52*np.ones(w.shape)
    
    # x, 31, 21
    
    arrs = (zero, first, second)
    narr = np.stack(arrs)
    d_arr = np.array([0, 100, 0])[:,None]
    R, T, A = optics.stack_calculation('s', narr, d_arr, 0.0, wreshape)
    
    fig, gs = wt.artists.create_figure(width='double', cols=[1,1])
    ax = plt.subplot(gs[0])
    ax.plot(w, R, label='R')
    ax.plot(w, A, label='A')
    ax.plot(w, A-A.min(), linewidth=2, label='shifted A')
    ax.plot(w, T, label='T')
    ax.set_xlim(1.3, 2.4)
    
    #ax.plot(d, channel='frac_abs', linewidth=2, label='shifted exp A')
    ax.grid()
    ax.legend()
    
    ax = plt.subplot(gs[1])
    ax.grid()
    ax.plot(w, out.imag, label='e imag')
    ax.plot(w, out.real, label='e real')
    ax.legend()
    
    