"""
Automated test code for eq

@author: Nick Gibbons
"""
import unittest
from numpy import array, zeros, absolute, log
from numpy.testing import assert_array_almost_equal
import eqc

Ru = 8.3144598 # uses gram-moles, or regular moles... to match M in lewis_library

class TestCEQ(unittest.TestCase):
    def test_pt(self):
        spnames = ['N2', 'O2', 'N', 'O', 'NO']
        T = 2500.0
        p = 0.1*101.35e3
        Xs0 = array([0.76, 0.23, 0.0, 0.0, 0.0])
        Xst= array([7.48543073e-01, 2.08366049e-01, 7.93468988e-07, 2.07645979e-02, 2.23254866e-02]) # Target with CEA mistake code
        #Xst = array([7.51474273e-01, 2.13161664e-01, 4.82204631e-07, 1.27384726e-02, 2.26251083e-02])


        eq = eqc.EqCalculator(spnames)

        Xs1 = eq.pt(p, T, Xs0, 0)
        Ys0 = eq.XtoY(Xs0)
        Ys1 = eq.XtoY(Xs1)
        Yst = eq.XtoY(Xst)
        assert_array_almost_equal(Ys1, Yst, decimal=6)
        #print("Done Y:   ", Ys1)
        #print("Target Y: ", Yst)
        return

    def test_rhou(self):
        spnames = ['N2', 'O2', 'N', 'O', 'NO']
        T = 2500.0
        p = 0.1*101.35e3
        Xs0 = array([0.76, 0.23, 0.0, 0.0, 0.0])
        Xst = array([7.48543073e-01, 2.08366049e-01, 7.93468988e-07, 2.07645979e-02, 2.23254866e-02])
        eq = eqc.EqCalculator(spnames)

        Mt = sum(Xst*eq.M)
        Rt = Ru/Mt
        rhot = p/Rt/T
        nt = 1/Mt
        cst = Xst*rhot/Mt
        nst = Xst/Mt # also cs/rhot 
        nt2 = nst.sum()

        ut = eq.get_u(Xst, T)
        ht = eq.get_h(Xst, T)

        Xs1, Teq = eq.rhou(rhot, ut, Xs0, 0)
        Ys0 = eq.XtoY(Xs0)
        Ys1 = eq.XtoY(Xs1)
        Yst = eq.XtoY(Xst)
        assert_array_almost_equal(Ys1, Yst, decimal=6)
        self.assertAlmostEqual(Teq, T, 4)
        return

    def test_ps(self):
        spnames = ['N2', 'O2', 'N', 'O', 'NO']
        T = 2500.0
        p = 0.1*101.35e3
        Xs0 = array([0.76, 0.23, 0.0, 0.0, 0.0])
        Xst = array([7.48543073e-01, 2.08366049e-01, 7.93468988e-07, 2.07645979e-02, 2.23254866e-02])
        eq = eqc.EqCalculator(spnames)

        Mt = sum(Xst*eq.M)
        Rt = Ru/Mt
        rhot = p/Rt/T
        nt = 1.0/Mt
        cst = Xst*rhot/Mt
        nst = Xst/Mt # also cs/rhot 
        nt2 = nst.sum()

        st = eq.get_s(Xst, T, p)
        pt = p

        Xs1, Teq = eq.ps(pt, st, Xs0, 0)
        Ys0 = eq.XtoY(Xs0)
        Ys1 = eq.XtoY(Xs1)
        Yst = eq.XtoY(Xst)
        assert_array_almost_equal(Ys1, Yst, decimal=6)
        self.assertAlmostEqual(Teq, T, 4)

if __name__=='__main__':
    unittest.main()
