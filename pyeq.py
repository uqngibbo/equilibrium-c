"""
Python wrapper for ceq core

@author: Nick Gibbons
"""

from string import ascii_letters
from numpy import array, zeros
from ctypes import cdll,c_double,POINTER,c_int
from lewis_thermo import get_species

letters = set(ascii_letters)
DBPATH = '/home/qungibbo/programs/us3d/1.0-RC22.12/props/lewis_thermo.db'

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

def pt(p, T, Xs0, nsp, nel, lewis, M, a):
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

    lib = cdll.LoadLibrary('./libceq.so')
    lib.pt(pp, Tp, Xs0p, nsp, nel, lewisp, Mp, ap, Xs1p)
    return Xs1

def rhou(rho, u, Xs0, nsp, nel, lewis, M, a):
    """ Call c library to compute equilibrium concentrations at fixed rho, u """
    Xs1 = zeros(Xs0.shape)
    rhop = c_double(rho)
    up = c_double(u)

    c_double_p = POINTER(c_double)
    Xs0p  = Xs0.ctypes.data_as(c_double_p)
    Mp    = M.ctypes.data_as(c_double_p)
    lewisp= lewis.ctypes.data_as(c_double_p)
    ap    = a.ctypes.data_as(c_double_p)
    Xs1p  = Xs1.ctypes.data_as(c_double_p)

    lib = cdll.LoadLibrary('./libceq.so')
    lib.rhou(rhop, up, Xs0p, nsp, nel, lewisp, Mp, ap, Xs1p)
    return Xs1


def test_pt():
    spnames = ['CO2', 'CO', 'O2']
    T = 2500.0
    p = 0.1*101.35e3
    Xs0 = array([1.0, 0.0, 0.0])
    Xst = array([0.66108962603325838,0.22594024931116111,0.11297012465558055])

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
    print("Computing")
    Xs1 = pt(p, T, Xs0, nsp, nel, lewisdata, M, a)
    print("Done: ", Xs1)
    print("Target: ", Xst)
    return

def test_rhou():
    spnames = ['CO2', 'CO', 'O2']
    T = 2500.0
    p = 0.1*101.35e3
    Xs0 = array([1.0, 0.0, 0.0])
    Xst = array([0.66108962603325838,0.22594024931116111,0.11297012465558055])

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

    Ru = 8.3144598 # uses gram-moles, or regular moles... to match M in lewis_library
    Mt = sum(Xst*M)
    Rt = Ru/Mt
    rhot = p/Rt/T
    nt = 1/Mt
    cst = Xst*rhot/Mt
    nst = Xst/Mt # also cs/rhot 
    nt2 = nst.sum()

    species = [get_species(sp) for sp in spnames]
    U0st= [(sp.H0onRT(T) - 1.0)*Ru*T for sp in species]
    ut = sum(njt*U0jt for njt,U0jt in zip(nst,U0st))

    print("a",a.flatten())
    print("ut",ut)
    print("nst",nst)
    print("rhot",rhot)
    print("Computing")
    Xs1 = rhou(rhot, ut, Xs0, nsp, nel, lewisdata, M, a)
    print("Done: ", Xs1)
    print("Target: ", Xst)
    return
    
if __name__=='__main__':
    #test_pt()
    test_rhou()
