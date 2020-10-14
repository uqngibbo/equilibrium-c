"""
Automated test code for ceq

@author: Nick Gibbons
"""
from numpy import array, zeros
import pyeq

def test_pt():
    spnames = 'N2,O2,C2H4,H2,CO,CO2,H2O'.split(',')
    T = 1024.278824
    p = 147514.463716
    Xs0= zeros(len(spnames))
    Xs0[0]=0.711202
    Xs0[1]=0.033864
    Xs0[2]=0.000046
    Xs0[3]=0.010316
    Xs0[4]=0.057133
    Xs0[5]=0.070316
    Xs0[6]=0.117124

    ceq = pyeq.EqCalculator(spnames)

    a = ceq.a
    M = ceq.M
    nel = ceq.nel
    nsp = ceq.nsp

    aij = a.reshape(nel,nsp).copy()
    Mmix = sum(Xs0*M)
    n0 = sum(Xs0/Mmix)
    ns0= Xs0/Mmix
    bi0 = (aij*ns0).sum(axis=1)
    print(spnames)
    print(ceq.elements)
    print("a\n", aij)
    print("bi0", bi0)

    print("Computing")
    Xs1 = ceq.pt(p, T, Xs0, 1)
    print("Done (Fixed conditions @ p={:1.6f} T={:1.6f})".format(p,T))
    print("Name  Init       Computed   Diff")
    for s,k in enumerate(spnames):
        print('{:>4}: {:1.6f}   {:1.6f}   {:1.6e}'.format(k, Xs0[s], Xs1[s], abs(Xs1[s]-Xs0[s])))

    Mmix = sum(Xs1*M)
    n1 = sum(Xs1/Mmix)
    ns1= Xs1/Mmix
    bi1 = (aij*ns1).sum(axis=1)
    print("bi0", bi0)
    print("bi1", bi1)
    return

if __name__=='__main__':
    test_pt()
    print(' ')
