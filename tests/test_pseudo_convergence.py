"""
Manual test code for eq

@author: Nick Gibbons
"""
import unittest
from numpy import array, zeros, absolute
from numpy.testing import assert_array_almost_equal
import eqc

class TestCEQ(unittest.TestCase):
    def test_pt(self, verbose=False):
        spnames = 'O2,CO,CO2'.split(',')
        T = 300.0
        p = 101e3
        Xs0= zeros(len(spnames))
        Xs0[0]=0.0
        Xs0[1]=0.0
        Xs0[2]=1.0
        Xst= zeros(len(spnames))
        Xst[0]=0.0
        Xst[1]=0.0
        Xst[2]=1.0

        eq = eqc.EqCalculator(spnames)
        M = eq.M
        a = eq.a
        nel = eq.nel
        nsp = eq.nsp
       
        aij = a.reshape(nel,nsp).copy()
        Mmix = sum(Xs0*M)
        n0 = sum(Xs0/Mmix)
        ns0= Xs0/Mmix
        bi0 = (aij*ns0).sum(axis=1)

        verbosity=0
        if verbose: verbosity=1
        Xs1 = eq.pt( p, T, Xs0, verbosity)
        if verbose:
            print("Name  Init       Computed   Diff")
            for s,k in enumerate(spnames):
                print('{:>4}: {:1.6f}   {:1.6f}   {:1.6e}'.format(k, Xs0[s], Xs1[s], abs(Xs1[s]-Xs0[s])))

        assert_array_almost_equal(Xs1, Xst, decimal=8)
        Mmix = sum(Xs1*M)
        n1 = sum(Xs1/Mmix)
        ns1= Xs1/Mmix
        bi1 = (aij*ns1).sum(axis=1)
        assert_array_almost_equal(bi0, bi1, decimal=8)
        return

if __name__=='__main__':
    unittest.main()
