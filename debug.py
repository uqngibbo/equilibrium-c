"""
Automated test code for ceq

@author: Nick Gibbons
"""

from numpy import array, zeros
from lewis_thermo import get_species
from pyeq import *

def test_pt():
    spnames = ['CO2', 'CO', 'O2']
    T = 500.0
    p = 0.1*101.35e3
    Xs0 = array([1.0, 0.0, 0.0])
    elements, nsp, nel, lewisdata, a, M = startup(spnames)

    print("Computing")
    Xs1 = pt(p, T, Xs0, nsp, nel, lewisdata, M, a, 1)
    print("Done: ", Xs1, Xs1.sum())
    return

if __name__=='__main__':
    test_pt()
