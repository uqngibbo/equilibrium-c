"""
Automated test code for eq

@author: Nick Gibbons
"""

import unittest
from numpy import array, zeros, absolute
from numpy.testing import assert_array_almost_equal
import eqc
from os import system, path
from re import search,split

def get_cea_test(p,T,spnames):
    #with open('ceatests/template.inp') as fp: text=fp.read()
    #text = text.replace('PRESSURE', str(p/101.35e3))
    #text = text.replace('TEMPERATURE', str(T))
    #with open('ceatests/co2.inp','w') as fp: fp.write(text)

    #system('cea2 ceatests/co2')

    filepath = path.dirname(__file__)
    cea_results = path.join(filepath, 'ceatests/co2.out')
    with open(cea_results) as fp: text=fp.read()

    regexstring = r'MOLE FRACTIONS\n\n( [\*]*[\S]+[\s]+[\s\.\d-]+\n)+\n'
    match = search(regexstring,text)
    if match==None: raise Exception("No match find in file!")

    lines = match.group().splitlines()
    keys = [i[:-8].strip().strip('*') for i in lines[2:-1]]
    values = [i[-8:].replace(' ','e').replace('-','e-') for i in lines[2:-1]]

    X = [0.0]*3
    for k,v in zip(keys,values):
        X[spnames.index(k)] = float(v)

    return X

class TestCEQ(unittest.TestCase):
    def test_pt(self, verbose=False):
        spnames = ['CO2', 'CO', 'O2']
        p = 0.1*101.35e3
        T = 2500.0
        Xs0 = array([1.0, 0.0, 0.0])

        eq = eqc.EqCalculator(spnames)
        Xcea1 = get_cea_test(p,T,spnames)
        Xs1 = eq.pt(p, T, Xs0)
        assert_array_almost_equal(Xs1, Xcea1, decimal=4)

if __name__=='__main__':
    unittest.main()
