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
    Ys0 = ceq.XtoY(Xs0)
    Ys1 = ceq.XtoY(Xs1)
    Yst = ceq.XtoY(Xst)
    print("Done:   ", Ys1)
    print("Target: ", Yst)
    return

def test_rhou(ceq, T, p, Xs0, Xst):
    Mt = sum(Xst*ceq.M)
    Rt = Ru/Mt
    rhot = p/Rt/T
    nt = 1/Mt
    cst = Xst*rhot/Mt
    nst = Xst/Mt # also cs/rhot 
    nt2 = nst.sum()

    ut = ceq.get_u(Xst, T)
    ht = ceq.get_h(Xst, T)

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

    st = ceq.get_s(Xst, T, p)
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
    spnames = ['N2', 'O2', 'N', 'O', 'NO']
    T = 2500.0
    p = 0.1*101.35e3
    Xs0 = array([0.76, 0.23, 0.0, 0.0, 0.0])
    Xst = array([0.74854377, 0.20836613, 0.0, 0.0207646, 0.0223255] )
    ceq = pyeq.EqCalculator(spnames)

    test_pt(ceq, T, p, Xs0, Xst)
    print(' ')
    test_rhou(ceq, T, p, Xs0, Xst)
    print(' ')
    test_ps(ceq, T, p, Xs0, Xst)
