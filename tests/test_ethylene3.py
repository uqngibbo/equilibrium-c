"""
Automated test code for eq

@author: Nick Gibbons
"""
import unittest
from numpy import array, zeros, absolute
from numpy.testing import assert_array_almost_equal
import eqc

class TestCEQ(unittest.TestCase):
    def test_pt(self, verbose=False):
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

        Xst= zeros(len(spnames))
        Xst[0]=7.3602e-1
        Xst[1]=0.0
        Xst[2]=0.0
        Xst[3]=0.0
        Xst[4]=0.0
        Xst[5]=1.3199e-1
        Xst[6]=1.3198e-1

        eq = eqc.EqCalculator(spnames)

        a = eq.a
        M = eq.M
        nel = eq.nel
        nsp = eq.nsp

        aij = a.reshape(nel,nsp).copy()
        Mmix = sum(Xs0*M)
        n0 = sum(Xs0/Mmix)
        ns0= Xs0/Mmix
        bi0 = (aij*ns0).sum(axis=1)

        verbosity=0
        if verbose: verbosity=1
        Xs1 = eq.pt(p, T, Xs0, verbosity)

        if verbose:
            print("Done (Fixed conditions @ p={:1.6f} T={:1.6f})".format(p,T))
            print("Name  Init       Computed   Diff")
            for s,k in enumerate(spnames):
                print('{:>4}: {:1.6f}   {:1.6f}   {:1.6e}'.format(k, Xs0[s], Xs1[s], abs(Xs1[s]-Xs0[s])))

        assert_array_almost_equal(Xs1, Xst, decimal=4)
        Mmix = sum(Xs1*M)
        n1 = sum(Xs1/Mmix)
        ns1= Xs1/Mmix
        bi1 = (aij*ns1).sum(axis=1)
        assert_array_almost_equal(bi0, bi1, decimal=8)
        return

if __name__=='__main__':
    unittest.main()
