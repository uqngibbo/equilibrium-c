"""
Automated test code for ceq

@author: Nick Gibbons
"""
from numpy import array, zeros
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

    #species = [get_species(sp) for sp in spnames]
    #U0st= [(sp.H0onRT(T) - 1.0)*Ru*T for sp in species]
    #ut = sum(njt*U0jt for njt,U0jt in zip(nst,U0st))
    ut = pyeq.get_u(lib, T, Xst, nsp, lewisdata, M)

    print("a",a.flatten())
    print("ut",ut)
    print("nst",nst)
    print("rhot",rhot)
    print("Computing")
    Xs1, Teq = pyeq.rhou(lib, rhot, ut, Xs0, nsp, nel, lewisdata, M, a, 1)
    print("Done X: ", Xs1)
    print("Target: ", Xst)
    print("Done T: ", T)
    print("Target: ", Teq)
    return

if __name__=='__main__':
    test_pt()
    print(' ')
    test_rhou()
