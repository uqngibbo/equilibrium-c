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
    print("Done Y:   ", Ys1)
    print("Target Y: ", Yst)
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
    print("Computing rhou....")
    Xs1, Teq = ceq.rhou(rhot, ut, Xs0, 1)
    Ys0 = ceq.XtoY(Xs0)
    Ys1 = ceq.XtoY(Xs1)
    Yst = ceq.XtoY(Xst)
    print("Done Y:   ", Ys1)
    print("Target Y: ", Yst)
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
    Ys0 = ceq.XtoY(Xs0)
    Ys1 = ceq.XtoY(Xs1)
    Yst = ceq.XtoY(Xst)
    print("Done Y:   ", Ys1)
    print("Target Y: ", Yst)
    print("Done T: ", Teq)
    print("Target: ", T)

if __name__=='__main__':
    spnames = ['N2', 'O2', 'N', 'O', 'NO', 'NO+', 'e-']
    T = 9000.0
    p = 0.1*101.35e3
    #Xs0 = array([0.75, 0.22, 0.0, 0.0, 0.0, 0.0, 0.0])
    Xs0 = array([1.2680e-3, 0.0, 7.7128e-1, 2.2490e-1, 2.8465e-5, 1.3e-3, 1.0e-3])
    Xs0/=(Xs0.sum())
    Xst = array([1.2680e-3, 0.0, 7.7128e-1, 2.2490e-1, 2.8465e-5, 1.2625e-3, 1.2625e-3])
    ceq = pyeq.EqCalculator(spnames)
    print(ceq.spnames)
    print(ceq.atoms)
    print(ceq.a)

    test_pt(ceq, T, p, Xs0, Xst)
    print(' ')
    #test_rhou(ceq, T, p, Xs0, Xst)
    #print(' ')
    #test_ps(ceq, T, p, Xs0, Xst)
