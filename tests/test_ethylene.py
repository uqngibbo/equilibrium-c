"""
Automated test code for eq

@author: Nick Gibbons
"""
import unittest
from numpy import array, zeros, absolute
from numpy.testing import assert_array_almost_equal
import eqc

def molef_from_massf(Y,Ms):
    Mmix = sum(Y/Ms)**-1
    return Y*Mmix/Ms

class TestCEQ(unittest.TestCase):
    def test_pt(self, verbose=False):
        spnames = 'N2,O2,C2H4,H2,CO,CO2,H2O'.split(',')
        eq = eqc.EqCalculator(spnames)
        Ms = eq.M
        ceaYs0 = array([0.75,  0.2, 0.05, 0.0, 0.0, 0.0, 0.0])
        ceaXst = {
            'CO': 1.2624e-2,
            'CO2': 8.9029e-2,
            'H2': 2.2749e-3,
            'H2O': 9.9379e-2,
            'N2': 7.6349e-1,
            'O2': 3.3208e-2
        }

        Xst = zeros(len(spnames))
        for k,v in ceaXst.items(): Xst[spnames.index(k)] = v

        Xs0 = molef_from_massf(ceaYs0, Ms)

        T = 2500.0
        p = 2.0*101.325e3

        verbosity = 0
        if verbose:
            print("Computing")
            verbosity=2
        Xs1 = eq.pt(p, T, Xs0, verbosity)

        if verbose:
            print("Done")
            print("Name  Init      Target    Computed")
            for s,k in enumerate(spnames):
                print('{:>4}: {:1.6f}  {:1.6f}  {:1.6f}'.format(k, Xs0[s], Xst[s], Xs1[s]))
        assert_array_almost_equal(Xs1, Xst, decimal=5)
        return

if __name__=='__main__':
    unittest.main()
