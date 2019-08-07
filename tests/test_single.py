"""
Automated test code for ceq

@author: Nick Gibbons
"""
from numpy import array, zeros, log
import pyeq

def test_pt():
    spnames = ['CO2', 'CO', 'O2']
    T = 2500.0
    p = 0.1*101.35e3
    Xs0 = array([1.0, 0.0, 0.0])
    Xst = array([0.66108962603325838,0.22594024931116111,0.11297012465558055])
    elements, nsp, nel, lewisdata, a, M = pyeq.startup(spnames)
    lib = pyeq.load_ceq_library()

    print("Computing")
    Xs1 = pyeq.pt(lib, p, T, Xs0, nsp, nel, lewisdata, M, a, 1)
    print("Done: ", Xs1)
    print("Target: ", Xst)
    return

def test_rhou():
    spnames = ['CO2', 'CO', 'O2']
    T = 2500.0
    p = 0.1*101.35e3
    Xs0 = array([1.0, 0.0, 0.0])
    Xst = array([0.66108962603325838,0.22594024931116111,0.11297012465558055])
    elements, nsp, nel, lewisdata, a, M = pyeq.startup(spnames)
    lib = pyeq.load_ceq_library()

    Ru = 8.3144598 # uses gram-moles, or regular moles... to match M in lewis_library
    Mt = sum(Xst*M)
    Rt = Ru/Mt
    rhot = p/Rt/T
    nt = 1/Mt
    cst = Xst*rhot/Mt
    nst = Xst/Mt # also cs/rhot 
    nt2 = nst.sum()

    ut = pyeq.get_u(lib, T, Xst, nsp, lewisdata, M)
    ht = pyeq.get_h(lib, T, Xst, nsp, lewisdata, M)

    print("a",a.flatten())
    print("ut",ut)
    print("ht",ht)
    print("nst",nst)
    print("rhot",rhot)
    print("Computing")
    Xs1, Teq = pyeq.rhou(lib, rhot, ut, Xs0, nsp, nel, lewisdata, M, a, 2)
    print("Done X: ", Xs1)
    print("Target: ", Xst)
    print("Done T: ", Teq)
    print("Target: ", T)
    return

def test_ps():
    spnames = ['CO2', 'CO', 'O2']
    T = 2500.0
    p = 0.1*101.35e3
    Xs0 = array([1.0, 0.0, 0.0])
    Xst = array([0.66108962603325838,0.22594024931116111,0.11297012465558055])
    elements, nsp, nel, lewisdata, a, M = pyeq.startup(spnames)
    lib = pyeq.load_ceq_library()

    Ru = 8.3144598 # uses gram-moles, or regular moles... to match M in lewis_library
    Mt = sum(Xst*M)
    Rt = Ru/Mt
    rhot = p/Rt/T
    nt = 1/Mt
    cst = Xst*rhot/Mt
    nst = Xst/Mt # also cs/rhot 
    nt2 = nst.sum()

    s0 = pyeq.get_s0(lib, T, Xst, nsp, lewisdata, M)
    st = s0 - (Ru*nst*(log(nst/nt) + log(p/1e5))).sum()
    pt = p

    print("a",a.flatten())
    print("st",st)
    print("nst",nst)
    print("pt",pt)
    print("Computing")
    Xs1, Teq = pyeq.ps(lib, pt, st, Xs0, nsp, nel, lewisdata, M, a, 1)
    print("Done X: ", Xs1)
    print("Target: ", Xst)
    print("Done T: ", Teq)
    print("Target: ", T)

if __name__=='__main__':
    test_pt()
    print(' ')
    test_rhou()
    print(' ')
    test_ps()
