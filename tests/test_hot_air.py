"""
Automated test code for eq

@author: Nick Gibbons
"""

import unittest
from numpy import array, zeros, absolute, linspace
import eqc
from gdtk.gas import GasModel, GasState

class TestCEQ(unittest.TestCase):
    def test_thermo(self, verbose=True):
        spnames = ['N2', 'O2', 'N', 'O', 'NO']
        T = 1900.0
        p = 100000.0
        Xs0 = array([0.76, 0.23, 0.0, 0.0, 0.0])
        eq = eqc.EqCalculator(spnames)
        Xs1 = eq.pt(p, T, Xs0)

        gmodel = GasModel('air-5sp-1T.lua')
        gs = GasState(gmodel)
        gs.p = p
        gs.T = T
        gs.massf = {k:Xsi for k,Xsi in zip(spnames, Xs1)}
        gs.update_thermo_from_pT()

        ceaoutput = """
        THERMODYNAMIC PROPERTIES
        
        P, BAR            1.0002
        T, K             1900.00
        RHO, KG/CU M    1.8324e-1
        H, KJ/KG         1858.81
        U, KJ/KG         1312.93
        G, KJ/KG        -15142.1
        S, KJ/(KG)(K)     8.9479
        
        M, (1/n)          28.940
        (dLV/dLP)t      -1.00004
        (dLV/dLT)p        1.0012
        Cp, KJ/(KG)(K)    1.3180
        GAMMAs            1.2795
        SON VEL,M/SEC      835.7
        
        MOLE FRACTIONS
        
        *NO             5.9082e-3
        *N2             7.6399e-1
        *O              1.4297e-4
        *O2             2.2996e-1
		*N              0.0000e00
        """
        ceadata = [i.strip().split('  ') for i in ceaoutput.splitlines()]
        ceadata = dict([(i[0], float(i[-1])) for i in ceadata if len(i)>2])
        ceaXs1 = array([ceadata['*N2'], ceadata['*O2'], ceadata['*N'], ceadata['*O'], ceadata['*NO']])
        eqdata = {}
        eqdata['h'] = eq.get_h(Xs1, T)
        eqdata['u'] = eq.get_u(Xs1, T)
        eqdata['s'] = eq.get_s(Xs1, T, p)
        eqdata['cp'] = eq.get_cp(Xs1, T)
        Teps = 1.0
        eqdata['cp2'] = (eq.get_h(Xs1, T+Teps) - eq.get_h(Xs1, T))/Teps

        Xs2 = eq.pt(p, T+Teps, Xs0)
        h2 = eq.get_h(Xs2, T+Teps)
        eqdata['cp3'] = (h2 - eqdata['h'])/Teps

        if verbose:
            print("          eqc    |   cea       |   gdtk ")
            print("XN2:    {:8f} |   {:8f}  | ".format(Xs1[0], ceaXs1[0]))
            print("XO2:    {:8f} |   {:8f}  |".format(Xs1[1], ceaXs1[1]))
            print("XN:     {:8f} |   {:8f}  |".format(Xs1[2], ceaXs1[2]))
            print("XO:     {:8f} |   {:8f}  |".format(Xs1[3], ceaXs1[3]))
            print("XNO:    {:8f} |   {:8f}  |".format(Xs1[4], ceaXs1[4]))
            print("h:      {:8.2f} |   {:8.2f}  | {:8.2f} (kJ/kg)".format(eqdata['h']/1000.0, ceadata['H, KJ/KG'], gmodel.enthalpy(gs)/1000.0))
            print("u:      {:8.2f} |   {:8.2f}  | {:8.2f} (kJ/kg)".format(eqdata['u']/1000.0, ceadata['U, KJ/KG'], gs.u/1000.0))
            print("s:      {:8f} |   {:8f}  | {:8f} (kJ/kg/K)".format(eqdata['s']/1000.0, ceadata['S, KJ/(KG)(K)'], gmodel.entropy(gs)/1000.0))
            print("cp:     {:8f} |   {:8f}  | {:8f}   |   {:8f}   | {:8f}   (kJ/kg/K)".format(
                   eqdata['cp']/1000.0, ceadata['Cp, KJ/(KG)(K)'], gs.Cp/1000.0, eqdata['cp2']/1000.0, eqdata['cp3']/1000.0))

        self.assertAlmostEqual(eqdata['h']/1000.0, ceadata['H, KJ/KG'], 0)
        self.assertAlmostEqual(eqdata['u']/1000.0, ceadata['U, KJ/KG'], 1)
        self.assertAlmostEqual(eqdata['s']/1000.0, ceadata['S, KJ/(KG)(K)'], 2)
        self.assertAlmostEqual(eqdata['cp3']/1000.0, ceadata['Cp, KJ/(KG)(K)'], 3)

if __name__=='__main__':
    unittest.main()
