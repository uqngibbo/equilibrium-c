"""
Automated test code for ceq

@author: Nick Gibbons
"""
import unittest
from numpy import array, zeros, absolute
from numpy.testing import assert_array_almost_equal
import pyeq


Ru = 8.3144598 # uses gram-moles, or regular moles... to match M in lewis_library

class TestCEQ(unittest.TestCase):
    def test_pt(self):
        spnames = ['N2', 'O2', 'N', 'O', 'NO', 'NO+', 'e-']
        T = 9000.0
        p = 0.1*101.35e3
        Xs0 = array([1.2680e-3, 0.0, 7.7128e-1, 2.2490e-1, 2.8465e-5, 1.3e-3, 1.0e-3])
        Xs0/=(Xs0.sum())
        Xst = array([1.2680e-3, 0.0, 7.7128e-1, 2.2490e-1, 2.8465e-5, 1.2625e-3, 1.2625e-3])
        ceq = pyeq.EqCalculator(spnames)

        Xs1 = ceq.pt(p, T, Xs0, 0)
        assert_array_almost_equal(Xs1, Xst, decimal=4)
        return

    def test_rhou(self):
        spnames = ['N2', 'O2', 'N', 'O', 'NO', 'NO+', 'e-']
        T = 9000.0
        p = 0.1*101.35e3
        Xs0 = array([1.2680e-3, 0.0, 7.7128e-1, 2.2490e-1, 2.8465e-5, 1.3e-3, 1.0e-3])
        Xs0/=(Xs0.sum())
        Xst = array([1.2680e-3, 0.0, 7.7128e-1, 2.2490e-1, 2.8465e-5, 1.2625e-3, 1.2625e-3])
        ceq = pyeq.EqCalculator(spnames)

        Mt = sum(Xst*ceq.M)
        Rt = Ru/Mt
        rhot = p/Rt/T
        nt = 1/Mt
        cst = Xst*rhot/Mt
        nst = Xst/Mt # also cs/rhot 
        nt2 = nst.sum()

        ut = ceq.get_u(Xst, T)

        Xs1, Teq = ceq.rhou(rhot, ut, Xs0, 0)
        assert_array_almost_equal(Xs1, Xst, decimal=4)
        self.assertAlmostEqual(Teq, T, 0)
        return

    def test_ps(self):
        spnames = ['N2', 'O2', 'N', 'O', 'NO', 'NO+', 'e-']
        T = 9000.0
        p = 0.1*101.35e3
        Xs0 = array([1.2680e-3, 0.0, 7.7128e-1, 2.2490e-1, 2.8465e-5, 1.3e-3, 1.0e-3])
        Xs0/=(Xs0.sum())
        Xst = array([1.2680e-3, 0.0, 7.7128e-1, 2.2490e-1, 2.8465e-5, 1.2625e-3, 1.2625e-3])
        ceq = pyeq.EqCalculator(spnames)

        Mt = sum(Xst*ceq.M)
        Rt = Ru/Mt
        rhot = p/Rt/T
        nt = 1.0/Mt
        cst = Xst*rhot/Mt
        nst = Xst/Mt # also cs/rhot 
        nt2 = nst.sum()

        s0 = ceq.get_s0(Xst, T)
        st = ceq.get_s(Xst, T, p)
        pt = p

        Xs1, Teq = ceq.ps(pt, st, Xs0, 0)
        assert_array_almost_equal(Xs1, Xst, decimal=4)
        self.assertAlmostEqual(Teq, T, 0)

if __name__=='__main__':
    unittest.main()
