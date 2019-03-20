"""
Python wrapper for ceq core

@author: Nick Gibbons
"""

from string import ascii_letters
from numpy import array, zeros
from ctypes import cdll,c_double,POINTER,c_int,byref

letters = set(ascii_letters)
#DBPATH = '/home/qungibbo/programs/us3d/1.0-RC22.12/props/lewis_thermo.db'
DBPATH = '/home/qungibbo/sourcecode/nasacea/thermo.inp'

def readdb(name):
    """ Retrieve species 'name' from the lewis_thermo.db file """
    with open(DBPATH) as fp:

        iterline=iter(fp) # Make an iterator for inner loop behaviour

        # Start looking for the beginning of a species entry
        for line in iterline:
            if not line[0] in letters: continue # Skip if not a letter in first position

            if line.startswith(name):        # This could be the one, check the rest
                header = line.strip().split()
                if header[0]!=name: continue # Nope false alarm

                lines = [line]               # We've found it!
                for nextline in iterline:
                    if nextline[0] in letters: break
                    lines.append(nextline)
                break                        # Break the outer for loop and jump to return lines
        else:
            raise Exception("Name: {} not found!".format(name))
    return lines

def get_thermo_data(name):
    """ Stripped down version of get_species from lewis_thermo.py """
    data = readdb(name)
    header = data[0].strip().split()
    if header[0]!=name:
        raise IOError("Database read failed! {}!={}".format(name, header[0]))

    info = data[1].strip().split()
    intervals = int(info[0])
    M = float(info[-2])/1e3 # kg/mol (regular mole, or gram-mole, unlike in cea_I)

    atomstrings = [data[1][10+8*i:10+8*(i+1)] for i in range(5)]
    atoms = {}
    for s in atomstrings:
        if s[0]==' ': continue
        name,amount = s.split()
        atoms[name] = float(amount)
    
    lewis = []
    for i in range(intervals):
        li = []
        line1 = data[3+3*i].replace('D','E') # Fix Fortran Double nonsense
        li.extend([float(line1[16*j:16*(j+1)]) for j in range(5)])

        line2 = data[4+3*i].replace('D','E') # Fix Fortran Double nonsense
        li.extend([float(line2[16*j:16*(j+1)]) for j in range(2)])
        li.extend([float(line2[16*j:16*(j+1)]) for j in range(3,5)])
        lewis.append(li)

    if len(lewis)!=3:
        lewis.append(lewis[-1])
    assert len(lewis)==3

    return atoms, M, lewis

def startup(spnames):
    atoms = []
    M = []
    lewisdata = []
    for sp in spnames:
        asp, Msp, lsp = get_thermo_data(sp)
        atoms.append(asp)
        M.append(Msp)
        lewisdata.append(lsp)

    M = array(M)
    lewisdata = array(lewisdata)

    elements = set()
    for a in atoms:
        for k in a.keys(): elements.add(k)
    elements = list(elements)
    elements.sort()
    print("elements", elements)

    nsp = len(spnames)
    nel = len(elements)

    a = zeros((len(elements),len(spnames)))
    for i,s in enumerate(atoms):
        for k,v in s.items():
            j = elements.index(k)
            a[j,i] = v

    print("a",a.flatten())
    return elements, nsp, nel, lewisdata, a, M

def load_ceq_library():
    """ Load the c library and set return types """
    lib = cdll.LoadLibrary('./libceq.so')
    lib.compute_Cp0_R.restype = c_double
    lib.pt.restype = c_int
    lib.rhou.restype = c_int
    lib.get_u.restype = c_double
    lib.batch_pt.restype = c_int
    lib.batch_rhou.restype = c_int
    lib.batch_u.restype = c_int
    return lib

def pt(lib, p, T, Xs0, nsp, nel, lewis, M, a,verbose=0):
    """ Call c library to compute equilibrium concentrations at fixed p, T """
    Xs1 = zeros(Xs0.shape)
    pp = c_double(p)
    Tp = c_double(T)

    c_double_p = POINTER(c_double)
    Xs0p  = Xs0.ctypes.data_as(c_double_p)
    Mp    = M.ctypes.data_as(c_double_p)
    lewisp= lewis.ctypes.data_as(c_double_p)
    ap    = a.ctypes.data_as(c_double_p)
    Xs1p  = Xs1.ctypes.data_as(c_double_p)

    recode = lib.pt(pp, Tp, Xs0p, nsp, nel, lewisp, Mp, ap, Xs1p, verbose)
    return Xs1

def rhou(lib, rho, u, Xs0, nsp, nel, lewis, M, a,verbose=0):
    """ Call c library to compute equilibrium concentrations at fixed rho, u """
    Xs1 = zeros(Xs0.shape)
    rhop = c_double(rho)
    up = c_double(u)
    Tp = c_double()
    Tref = byref(Tp)

    c_double_p = POINTER(c_double)
    Xs0p  = Xs0.ctypes.data_as(c_double_p)
    Mp    = M.ctypes.data_as(c_double_p)
    lewisp= lewis.ctypes.data_as(c_double_p)
    ap    = a.ctypes.data_as(c_double_p)
    Xs1p  = Xs1.ctypes.data_as(c_double_p)

    recode = lib.rhou(rhop, up, Xs0p, nsp, nel, lewisp, Mp, ap, Xs1p, Tref, verbose)
    T = Tp.value
    return Xs1, T

def get_u(lib, T, X, nsp, lewis, M):
    """ Call c library to compute internal energy at fixed composition and temperature """
    Tp = c_double(T)

    c_double_p = POINTER(c_double)
    Xp    = X.ctypes.data_as(c_double_p)
    Mp    = M.ctypes.data_as(c_double_p)
    lewisp= lewis.ctypes.data_as(c_double_p)

    u = lib.get_u(Tp, Xp, nsp, lewisp, Mp)
    return u

def batch_pt(lib, p, T, Xs0, nsp, nel, lewis, M, a,verbose=0):
    """ Call c library to compute equilibrium concentrations at fixed p, T """
    N, nspcheck = Xs0.shape
    if nspcheck!=nsp: raise Exception("nsp ({}) != Xs0.shape[1] ({})".format(nsp, nspcheck))
    if N!=p.size: raise Exception("p.size ({}) != Xs0.shape[0] ({})".format(p.size, N))
    if N!=T.size: raise Exception("T.size ({}) != Xs0.shape[0] ({})".format(T.size, N))

    Xs1 = zeros(Xs0.shape)

    c_double_p = POINTER(c_double)
    pp    = p.ctypes.data_as(c_double_p)
    Tp    = T.ctypes.data_as(c_double_p)
    Xs0p  = Xs0.ctypes.data_as(c_double_p)
    Mp    = M.ctypes.data_as(c_double_p)
    lewisp= lewis.ctypes.data_as(c_double_p)
    ap    = a.ctypes.data_as(c_double_p)
    Xs1p  = Xs1.ctypes.data_as(c_double_p)

    recode = lib.batch_pt(N, pp, Tp, Xs0p, nsp, nel, lewisp, Mp, ap, Xs1p, verbose)
    return Xs1

def batch_rhou(lib, rho, u, Xs0, nsp, nel, lewis, M, a,verbose=0):
    """ Call c library to compute equilibrium concentrations at fixed rho, u """
    N, nspcheck = Xs0.shape
    if nspcheck!=nsp: raise Exception("nsp ({}) != Xs0.shape[1] ({})".format(nsp, nspcheck))
    if N!=rho.size: raise Exception("rho.size ({}) != Xs0.shape[0] ({})".format(rho.size, N))
    if N!=u.size: raise Exception("u.size ({}) != Xs0.shape[0] ({})".format(u.size, N))

    Xs1 = zeros(Xs0.shape)
    T = zeros(rho.shape)

    c_double_p = POINTER(c_double)
    rhop  = rho.ctypes.data_as(c_double_p)
    up    = u.ctypes.data_as(c_double_p)
    Tp    = T.ctypes.data_as(c_double_p)
    Xs0p  = Xs0.ctypes.data_as(c_double_p)
    Mp    = M.ctypes.data_as(c_double_p)
    lewisp= lewis.ctypes.data_as(c_double_p)
    ap    = a.ctypes.data_as(c_double_p)
    Xs1p  = Xs1.ctypes.data_as(c_double_p)

    recode = lib.batch_rhou(N, rhop, up, Xs0p, nsp, nel, lewisp, Mp, ap, Xs1p, Tp, verbose)
    return Xs1, T

if __name__=='__main__':
    print("Called pyeq main!")
