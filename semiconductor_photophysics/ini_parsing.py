import WrightTools as wt

def get_BK_numerical_params(ini, section_name='numerical_params'):
    xmax = ini.read(section_name, 'xmax')
    xnum = ini.read(section_name, 'xnum')
    nmax = ini.read(section_name, 'nmax')
    return (xmax, xnum, nmax)

def get_BKs(ini):
    out = []
    sections = ini.sections
    for section in sections:
        if 'BK' in section:
            tup = get_BK(ini, section)
            out.append(tup)
    return out

def get_BK(ini, section):
    Eg0    = ini.read(section, 'Eg0')
    G      = ini.read(section, 'G')
    a0     = ini.read(section, 'a0')
    k      = ini.read(section, 'k')
    T      = ini.read(section, 'T')
    rcv    = ini.read(section, 'rcv')
    mu_e   = ini.read(section, 'mu_e')
    mu_h   = ini.read(section, 'mu_h')
    m_star = ini.read(section, 'm_star')
    return (Eg0, G, a0, k, T, rcv, mu_e, mu_h, m_star)
    
def get_Lors(ini):
    out = []
    sections = ini.sections
    for section in sections:
        if 'Lorentzian' in section:
            tup = get_lor(ini, section)
            out.append(tup)
    return out

def get_lor(ini, section):
    E0 = ini.read(section, 'E0')
    G  = ini.read(section, 'G')
    A  = ini.read(section, 'A')
    A0 = ini.read(section, 'A0')
    return (E0, G, A, A0)

def read_full_sim_params(p):
    ini = wt.kit.INI(p)
    num_params = get_BK_numerical_params(ini, section_name='numerical_params')
    BKs = get_BKs(ini)
    Lors = get_Lors(ini)
    out = {}
    out['num_params'] = num_params
    out['BKs'] = BKs
    out['Lors'] = Lors
    return out
    
