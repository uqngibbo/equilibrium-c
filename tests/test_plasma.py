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
        spnames = ['N2', 'O2', 'N', 'O', 'NO', 'NO+', 'e-']
        T = 9000.0
        p = 0.1*101.35e3
        Xs0 = array([1.2680e-3, 0.0, 7.7128e-1, 2.2490e-1, 2.8465e-5, 1.3e-3, 1.0e-3])
        Xs0/=(Xs0.sum())
        Xst= array([1.2680e-3, 0.0,       7.7128e-1, 2.2490e-1, 2.8465e-5, 1.2625e-3, 1.2625e-3])
        #Xst = array([3.4238e-3, 1.1395e-6, 7.6860e-1, 2.2536e-1, 7.7288e-5, 1.2616e-3, 1.2616e-3])

        eq = eqc.EqCalculator(spnames)

        Xs1 = eq.pt(p, T, Xs0, 0)
        assert_array_almost_equal(Xs1, Xst, decimal=4)
        return

    def test_rhou(self):
        spnames = ['N2', 'O2', 'N', 'O', 'NO', 'NO+', 'e-']
        T = 9000.0
        p = 0.1*101.35e3
        Xs0 = array([1.2680e-3, 0.0, 7.7128e-1, 2.2490e-1, 2.8465e-5, 1.3e-3, 1.0e-3])
        Xs0/=(Xs0.sum())
        Xst = array([1.2680e-3, 0.0, 7.7128e-1, 2.2490e-1, 2.8465e-5, 1.2625e-3, 1.2625e-3])
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
        assert_array_almost_equal(Xs1, Xst, decimal=4)
        self.assertAlmostEqual(Teq, T, 0)

    def test_11_sp(self):
        spnames = ['N2', 'N2+', 'NO', 'NO+', 'O2', 'O2+', 'N', 'N+', 'O', 'O+', 'e-']
        T = 9000.0
        p = 0.1*101.35e3
        Xs0 = array([0.77, 0.0, 0.0, 0.0, 0.23, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        Xs0/=(Xs0.sum())
        eq = eqc.EqCalculator(spnames)
        Xst = array([1.12808714e-03, 1.66990829e-05, 2.62323131e-05, 5.67521824e-05,
                     3.98415990e-07, 1.05942060e-07, 7.27382202e-01, 2.12663187e-02,
                     2.19710628e-01, 4.53634979e-03, 2.58762257e-02])
        #Xst = array([3.18066841e-03, 2.82221200e-05, 7.40990843e-05, 9.60904185e-05,
        #             1.12749213e-06, 1.79707597e-07, 7.40803974e-01, 1.29823770e-02,
        #             2.24177594e-01, 2.77439912e-03, 1.58812684e-02])


        Xs1 = eq.pt(p, T, Xs0, 0)
        #print("Ys1: ", eq.XtoY(Xs1))
        #print("Xs1: ", dict(zip(spnames, Xs1)))
        assert_array_almost_equal(Xs1, Xst, decimal=8)

if __name__=='__main__':
    unittest.main()
