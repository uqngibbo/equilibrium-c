"""
Automated test code for eq

@author: Nick Gibbons
"""
import unittest
from numpy import array, zeros, absolute
from numpy.testing import assert_array_almost_equal
import eqc

Ru = 8.3144598 # uses gram-moles, or regular moles... to match M in lewis_library

class TestCEQ(unittest.TestCase):
    def test_pt(self):
        spnames = ['CO2', 'CO', 'O2']
        T = 2500.0
        p = 0.1*101.35e3
        Xs0 = array([1.0, 0.0, 0.0])
        #Xst = array([0.66108962603325838,0.22594024931116111,0.11297012465558055])
        Xst = array([0.66010397,  0.22659735,  0.11329868]) # Old code with CEA mistake

        eq = eqc.EqCalculator(spnames)

        Xs1 = eq.pt(p, T, Xs0, 0)
        assert_array_almost_equal(Xs1, Xst, decimal=8)
        return

    def test_rhou(self):
        spnames = ['CO2', 'CO', 'O2']
        T = 2500.0
        p = 0.1*101.35e3
        Xs0 = array([1.0, 0.0, 0.0])
        Xst = array([0.66010397,  0.22659735,  0.11329868])
        eq = eqc.EqCalculator(spnames)

        Mt = sum(Xst*eq.M)
        Rt = Ru/Mt
        rhot = p/Rt/T
        nt = 1/Mt
        cst = Xst*rhot/Mt
        nst = Xst/Mt # also cs/rhot 
        nt2 = nst.sum()

        ut = eq.get_u(Xst, T)

        Xs1, Teq = eq.rhou(rhot, ut, Xs0, 0)
        assert_array_almost_equal(Xs1, Xst, decimal=8)
        self.assertAlmostEqual(Teq, T, 4)
        return

    def test_ps(self):
        spnames = ['CO2', 'CO', 'O2']
        T = 2500.0
        p = 0.1*101.35e3
        Xs0 = array([1.0, 0.0, 0.0])
        Xst = array([0.66010397,  0.22659735,  0.11329868])
        eq = eqc.EqCalculator(spnames)

        Mt = sum(Xst*eq.M)
        Rt = Ru/Mt
        rhot = p/Rt/T
        nt = 1.0/Mt
        cst = Xst*rhot/Mt
        nst = Xst/Mt # also cs/rhot 
        nt2 = nst.sum()

        s0 = eq.get_s0(Xst, T)
        st = eq.get_s(Xst, T, p)
        pt = p

        Xs1, Teq = eq.ps(pt, st, Xs0, 0)
        assert_array_almost_equal(Xs1, Xst, decimal=8)
        self.assertAlmostEqual(Teq, T, 4)

    def test_rhot(self):
        spnames = ['CO2', 'CO', 'O2']
        T = 2500.0
        p = 0.1*101.35e3
        Xs0 = array([1.0, 0.0, 0.0])
        Xst = array([0.66010397,  0.22659735,  0.11329868])
        eq = eqc.EqCalculator(spnames)

        Mt = sum(Xst*eq.M)
        Rt = Ru/Mt
        rhot = p/Rt/T

        Xs1 = eq.rhot(rhot, T, Xs0, 0)
        assert_array_almost_equal(Xs1, Xst, decimal=7) # Seven?
        return

if __name__=='__main__':
    unittest.main()
