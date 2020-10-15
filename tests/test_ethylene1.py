"""
Test code for ceq: near zero trace species

@author: Nick Gibbons
"""

import unittest
from numpy import array, zeros, absolute
from numpy.testing import assert_array_almost_equal
import pyeq

class TestCEQ(unittest.TestCase):
    def test_pt(self, verbose=False):
        spnames = 'N2,O2,C2H4,H2,CO,CO2,H2O'.split(',')

        T = 300.0
        p = 100000.0
        Xs0= zeros(len(spnames))
        Xs0[0]=0.79
        Xs0[1]=0.21
        Xs0[2]=0.0
        Xs0[3]=0.0
        Xs0[4]=0.0
        Xs0[5]=0.0
        Xs0[6]=0.0

        Xst= zeros(len(spnames))
        Xst[0]=0.79
        Xst[1]=0.21
        Xst[2]=0.0
        Xst[3]=0.0
        Xst[4]=0.0
        Xst[5]=0.0
        Xst[6]=0.0
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

        verbosity = 0
        if verbose: verbosity=2
        Xs1 = ceq.pt(p, T, Xs0, verbosity)
        
        if verbose:
            print("Name  Init       Target     Computed   Diff")
            for s,k in enumerate(spnames):
                print('{:>4}: {:1.6f}   {:1.6f}   {:1.6f}   {:1.6e}'.format(k, Xs0[s], Xst[s], Xs1[s], abs(Xs1[s]-Xst[s])))

        assert_array_almost_equal(Xs1, Xst, decimal=5)
        Mmix = sum(Xs1*M)
        n1 = sum(Xs1/Mmix)
        ns1= Xs1/Mmix
        bi1 = (aij*ns1).sum(axis=1)
        assert_array_almost_equal(bi0, bi1, decimal=8)
        return

if __name__=='__main__':
    unittest.main()
