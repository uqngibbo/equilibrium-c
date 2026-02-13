"""
Test the saturated pt capability by computing Bprime at an ablating surface.

@author: Nick Gibbons
"""

import unittest
from numpy import array, linspace, zeros
from numpy.testing import assert_array_almost_equal
import eqc

def calc_Bprime(Ys, spnames, eq):
    # We need to know the mass fraction of each element in the mixture.
    Ns = {
        'N':zeros(Ys['N2'].shape),
        'O':zeros(Ys['N2'].shape),
        'C':zeros(Ys['N2'].shape),
    }

    for s,name in enumerate(spnames):
        for j,el in enumerate(eq.elements):
            Ns[el] += eq.a[j,s]*Ys[name]*eq.M[spnames.index(el)]/eq.M[s]

    Le = 1.0
    Bprime = Ns['C']/(Ns['N']+Ns['O'])/Le
    return Bprime

class TestCEQ(unittest.TestCase):
    def run_sat_calc(self, p, T):
        eq_carbon = eqc.EqCalculator(['C(gr)'])
        Xc = array([1.0])
        gc0 = eq_carbon.get_h(Xc, T) - T*eq_carbon.get_s0(Xc, T)
        Gc0 = gc0*eq_carbon.M[0]

        spnames = ['N2', 'O2', 'CO2', 'CO', 'N', 'O',
                   'C2', 'C3', 'C4', 'C5', 'CN', 'CNC',
                   'C2N2', 'C4N2', 'C']
        eq = eqc.EqCalculator(spnames)
        ic = eq.elements.index('C')

        bC = 0.7
        bO = 1.5

        # See derivation from 6/02/26, assumes nCO=nCO2 because reasons
        nO2  = (bO - 1.5*bC)/2.0
        nCO2 = (bO - 2.0*nO2)/3.0
        nCO = nCO2
        nN2 = 0.767/0.233*(2*nO2 + 2*nCO2 + 1*nCO)/2.0

        nsp = len(spnames)

        Xs0_init = zeros(nsp)
        Xs0_init[spnames.index('N2')]  = nN2
        Xs0_init[spnames.index('O2')]  = nO2
        Xs0_init[spnames.index('CO2')] = nCO2
        Xs0_init[spnames.index('CO')]  = nCO

        Xs0 = Xs0_init/Xs0_init.sum()
        Xs1 = eq.satpt(p, T, Gc0, ic, Xs0, 0)
        Ys1 = eq.XtoY(Xs1)
        Ys = {}
        for isp,sp in enumerate(spnames):
            Ys[sp] = Ys1[isp]

        Bprime = calc_Bprime(Ys, spnames, eq)
        return Bprime

    def test_satpt_CO2(self):
        """ At lowish temps, CO2 is in equilibrium with the surface """
        p = 101.35e3
        T = 600.0
        Bprime = self.run_sat_calc(p, T)
        Bprime_target = 0.0968242
        self.assertAlmostEqual(Bprime, Bprime_target, 4)

    def test_satpt_CO(self):
        """ At medium temps, we get another plateau associated with CO """
        p = 101.35e3
        T = 2000.0
        Bprime = self.run_sat_calc(p, T)
        Bprime_target = 0.1933858
        self.assertAlmostEqual(Bprime, Bprime_target, 4)

    def test_satpt_C(self):
        """ At high temps, the B' rises rapidly as the carbon begins to sublime """
        p = 101.35e3
        T = 3800.0
        Bprime = self.run_sat_calc(p, T)
        Bprime_target = 1.022456
        self.assertAlmostEqual(Bprime, Bprime_target, 4)

if __name__=='__main__':
    unittest.main()
