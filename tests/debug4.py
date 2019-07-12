"""
Test code for ceq: near zero trace species

@author: Nick Gibbons
"""
from numpy import array, zeros
import pyeq

def test_pt():
    spnames = 'N2,O2,C2H4,H2,CO,CO2,H2O'.split(',')
    T = 1253.579055
    p = 154800.574538
    Xs0= zeros(len(spnames))
    Xs0[0]=7.882249e-01
    Xs0[1]=2.114530e-01
    Xs0[2]=1.683690e-10
    Xs0[3]=5.046495e-12
    Xs0[4]=5.081258e-12
    Xs0[5]=3.220202e-04
    Xs0[6]=4.257038e-14

    #Xs0[0]=7.882249e-01
    #Xs0[1]=2.114530e-01
    #Xs0[2]=0.0
    #Xs0[3]=0.0
    #Xs0[4]=0.0
    #Xs0[5]=3.220202e-04
    #Xs0[6]=0.0


    elements, nsp, nel, lewisdata, a, M = pyeq.startup(spnames)

    lib = pyeq.load_ceq_library()

    aij = a.reshape(nel,nsp).copy()
    Mmix = sum(Xs0*M)
    n0 = sum(Xs0/Mmix)
    ns0= Xs0/Mmix
    bi0 = (aij*ns0).sum(axis=1)
    print(spnames)
    print(elements)
    print("a\n", aij)
    print("bi0", bi0)

    print("Computing")
    Xs1 = pyeq.pt(lib, p, T, Xs0, nsp, nel, lewisdata, M, a, 1)
    print("Done")
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
