"""
Automated test code for ceq

@author: Nick Gibbons
"""
from numpy import array, zeros, log
import pyeq

Ru = 8.3144598 # uses gram-moles, or regular moles... to match M in lewis_library

def test_pt(ceq, T, p, Xs0, Xst):
    print("Computing")
    Xs1 = ceq.pt(p, T, Xs0, 1)
    print("Done: ", Xs1)
    print("Target: ", Xst)
    return

def test_rhou(ceq, T, p, Xs0, Xst):
    Mt = sum(Xst*ceq.M)
    Rt = Ru/Mt
    rhot = p/Rt/T
    nt = 1/Mt
    cst = Xst*rhot/Mt
    nst = Xst/Mt # also cs/rhot 
    nt2 = nst.sum()

    ut = ceq.get_u(T, Xst)
    ht = ceq.get_h(T, Xst)

    print("ut",ut)
    print("ht",ht)
    print("nst",nst)
    print("rhot",rhot)
    print("Computing")
    Xs1, Teq = ceq.rhou(rhot, ut, Xs0, 1)
    print("Done X: ", Xs1)
    print("Target: ", Xst)
    print("Done T: ", Teq)
    print("Target: ", T)
    return

def test_ps(ceq, T, p, Xs0, Xst):
    Mt = sum(Xst*ceq.M)
    Rt = Ru/Mt
    rhot = p/Rt/T
    nt = 1.0/Mt
    cst = Xst*rhot/Mt
    nst = Xst/Mt # also cs/rhot 
    nt2 = nst.sum()

    s0 = ceq.get_s0(T, Xst)
    st = s0 - (Ru*nst*(log(nst/nt) + log(p/1e5))).sum()
    pt = p

    print("st",st)
    print("nst",nst)
    print("pt",pt)
    print("Computing")
    Xs1, Teq = ceq.ps(pt, st, Xs0, 1)
    print("Done X: ", Xs1)
    print("Target: ", Xst)
    print("Done T: ", Teq)
    print("Target: ", T)

if __name__=='__main__':
    spnames = ['CO2', 'CO', 'O2']
    T = 2500.0
    p = 0.1*101.35e3
    Xs0 = array([1.0, 0.0, 0.0])
    Xst = array([0.66108962603325838,0.22594024931116111,0.11297012465558055])
    ceq = pyeq.EqCalculator(spnames)

    test_pt(ceq, T, p, Xs0, Xst)
    print(' ')
    test_rhou(ceq, T, p, Xs0, Xst)
    print(' ')
    test_ps(ceq, T, p, Xs0, Xst)
