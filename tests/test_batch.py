"""
Automated test code for ceq

@author: Nick Gibbons
"""
import unittest
from numpy import array, zeros, absolute
from numpy.testing import assert_array_almost_equal
import pyeq

class TestCEQ(unittest.TestCase):
    def test_pt(self):
        spnames = ['CO2', 'CO', 'O2']
        ceq = pyeq.EqCalculator(spnames)

        T = array([2500.0, 2000.0])
        p = array([0.1*101.35e3, 1.0*101.35e3])
        Xs0 = array([[1.0, 0.0, 0.0],
                     [1.0, 0.0, 0.0]])
        Xst = zeros(Xs0.shape)
        np,ns = Xs0.shape

        for i in range(np):
            Xs0i = Xs0[i].copy()
            Ti = T[i]
            pi = p[i]
            Xs1i = ceq.pt(pi, Ti, Xs0i)
            Xst[i,:] = Xs1i[:]

        Xs1 = ceq.batch_pt(p, T, Xs0, 0)
        assert_array_almost_equal(Xs1, Xst)
        return

    def test_rhou(self):
        spnames = ['CO2', 'CO', 'O2']
        ceq = pyeq.EqCalculator(spnames)

        u  = array([-5e6, -4.5e6])
        rho= array([0.019, 0.025])
        Xs0 = array([[1.0, 0.0, 0.0],
                     [1.0, 0.0, 0.0]])
        Xst = zeros(Xs0.shape)
        Tt = zeros(rho.shape)
        np,ns = Xs0.shape

        for i in range(np):
            Xs0i = Xs0[i].copy()
            ui = u[i]
            rhoi = rho[i]
            Xs1i, Ti = ceq.rhou(rhoi, ui, Xs0i)
            Xst[i,:] = Xs1i[:]
            Tt[i] = Ti

        Xs1, T = ceq.batch_rhou(rho, u, Xs0, 0)
        assert_array_almost_equal(Xs1, Xst)
        assert_array_almost_equal(T, Tt)
        return

    # FIXME: We don't have this hooked up yet.
    #def test_ps(self):
    #    spnames = ['CO2', 'CO', 'O2']
    #    ceq = pyeq.EqCalculator(spnames)

    #    p = array([10135.0, 12000.0])
    #    s = array([8484.85, 9000.0])
    #    Xs0 = array([[1.0, 0.0, 0.0],
    #                 [1.0, 0.0, 0.0]])
    #    Xst = zeros(Xs0.shape)
    #    Tt = zeros(s.shape)
    #    np,ns = Xs0.shape

    #    for i in range(np):
    #        Xs0i = Xs0[i].copy()
    #        pi = p[i]
    #        si = s[i]
    #        Xs1i, Ti = ceq.ps(pi, si, Xs0i)
    #        Xst[i,:] = Xs1i[:]
    #        Tt[i] = Ti

    #    Xs1, T = ceq.batch_ps(p, s, Xs0, 0)
    #    assert_array_almost_equal(Xs1, Xst)
    #    assert_array_almost_equal(T, Tt)
    #    return

    def test_u(self):
        spnames = ['CO2', 'CO', 'O2']
        ceq = pyeq.EqCalculator(spnames)

        T = array([2500.0, 2200.0])
        Xs0 = array([[0.7, 0.2, 0.1],
                     [0.7, 0.2, 0.1]])
        ut = zeros(T.shape)
        np,ns = Xs0.shape

        for i in range(np):
            Xs0i = Xs0[i].copy()
            Ti = T[i]
            ut[i] = ceq.get_u(Xs0i, Ti)

        u = ceq.batch_u(Xs0, T)
        assert_array_almost_equal(u, ut)
        return

if __name__=='__main__':
    unittest.main()
