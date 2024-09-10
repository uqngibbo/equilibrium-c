"""
Test for the Lagrangian derivatives being zero.

@author: Nick Gibbons
"""
import unittest
from numpy import array, zeros, absolute, log, sqrt
import eqc

class TestCEQ(unittest.TestCase):
    def test_pt_lagrangian(self):
        spnames = ['N2', 'O2', 'N', 'O', 'NO']
        T = 2500.0
        p = 0.1*101.35e3
        Xs0 = array([0.767, 0.233, 0.0, 0.0, 0.0])

        eq = eqc.EqCalculator(spnames)
        Xs1 = eq.pt(p, T, Xs0, 0)

        dLdn = eq.verify_equilibrium(p, T, Xs1, 0)
        dLdn_L2 = sqrt(sum(i*i for i in dLdn))

        # The numerical derivatives are quite noisy,
        # especially considering L is about 1e7. So 
        # within one decimal place of 0.0 is a pretty
        # good test.
        self.assertAlmostEqual(dLdn_L2, 0.0, 1)

if __name__=='__main__':
    unittest.main()
