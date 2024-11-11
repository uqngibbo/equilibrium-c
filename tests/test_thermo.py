"""
Automated test code for eq

@author: Nick Gibbons
"""

import unittest
from numpy import array, zeros, absolute, linspace
import eqc

class TestCEQ(unittest.TestCase):
    def test_thermo(self, verbose=False):
        p = 1.0132*1e5
        T = 3698.04
        spnames = ['CO2', 'CO', 'O2']
        Xs0 = array([1.0, 0.0, 0.0])
        eq = eqc.EqCalculator(spnames)
        Xs1 = eq.pt(p, T, Xs0)

        ceaoutput = """
        THERMODYNAMIC PROPERTIES

        P, BAR            1.0132
        T, K             3698.04
        RHO, KG/CU M   1.0265e-1
        H, KJ/KG          576.89
        U, KJ/KG         -410.21
        G, KJ/KG        -35107.3
        S, KJ/(KG)(K)     9.6495



        M, (1/n)          31.149
        (dLV/dLP)t      -1.02399
        (dLV/dLT)p        1.4203
        Cp, KJ/(KG)(K)    3.3065
        GAMMAs            1.1612
        SON VEL,M/SEC     1070.6

        MOLE FRACTIONS

        *CO            5.8443e-1
        *CO2           1.2336e-1
        *O2            2.9221e-1
        """
        ceadata = [i.strip().split('  ') for i in ceaoutput.splitlines()]
        ceadata = dict([(i[0], float(i[-1])) for i in ceadata if len(i)>2])
        ceaXs1 = array([ceadata['*CO2'], ceadata['*CO'], ceadata['*O2']])
        eqdata = {}
        eqdata['h'] = eq.get_h(Xs1, T)
        eqdata['u'] = eq.get_u(Xs1, T)
        eqdata['s'] = eq.get_s(Xs1, T, p)

        if verbose:
            print("          eq    |   cea")
            print("XCO2:   {:8f} |   {:8f} ".format(Xs1[0], ceaXs1[0]))
            print("XCO:    {:8f} |   {:8f} ".format(Xs1[1], ceaXs1[1]))
            print("XO2:    {:8f} |   {:8f} ".format(Xs1[2], ceaXs1[2]))
            print("h:    {:8f} | {:8f}  (kJ/kg)".format(eqdata['h']/1000.0, ceadata['H, KJ/KG']))
            print("u:   {:8f} |{:8f}  (kJ/kg)".format(eqdata['u']/1000.0, ceadata['U, KJ/KG']))
            print("s:      {:8f} |   {:8f}  (kJ/kg/K)".format(eqdata['s']/1000.0, ceadata['S, KJ/(KG)(K)']))

        self.assertAlmostEqual(eqdata['h']/1000.0, ceadata['H, KJ/KG'], 1)
        self.assertAlmostEqual(eqdata['u']/1000.0, ceadata['U, KJ/KG'], 1)
        self.assertAlmostEqual(eqdata['s']/1000.0, ceadata['S, KJ/(KG)(K)'], 4)

if __name__=='__main__':
    unittest.main()
