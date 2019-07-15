"""
Automated test code for ceq

@author: Nick Gibbons
"""
from numpy import array, zeros
import pyeq

def test_pt():
    spnames = ['CO2', 'CO', 'O2']
    T = array([2500.0, 2000.0])
    p = array([0.1*101.35e3, 1.0*101.35e3])
    Xs0 = array([[1.0, 0.0, 0.0],
                 [1.0, 0.0, 0.0]])
    Xst = array([0.66108962603325838,0.22594024931116111,0.11297012465558055])
    elements, nsp, nel, lewisdata, a, M = pyeq.startup(spnames)
    lib = pyeq.load_ceq_library()

    print("Computing")
    Xs1 = pyeq.batch_pt(lib, p, T, Xs0, nsp, nel, lewisdata, M, a, 1)
    print("Done: ", Xs1)
    print("Target: ", Xst)
    return

def test_u():
    spnames = ['CO2', 'CO', 'O2']
    T = array([2500.0, 2200.0])
    Xs0 = array([[0.7, 0.2, 0.1],
                 [0.7, 0.2, 0.1]])
    elements, nsp, nel, lewisdata, a, M = pyeq.startup(spnames)
    lib = pyeq.load_ceq_library()

    print("Computing")
    u = pyeq.batch_u(lib, T, Xs0, nsp, lewisdata, M)
    print("Done: ", u)
    return

if __name__=='__main__':
    test_pt()
    print(' ')
    test_u()
