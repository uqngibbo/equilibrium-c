"""
Automated test code for ceq

@author: Nick Gibbons
"""
from numpy import array, zeros
import pyeq

def test_pt(ceq, p, T, X0, Xst):

    print("Computing")
    Xs1 = ceq.batch_pt(p, T, Xs0, 1)
    print("Done: ", Xs1)
    print("Target: ", Xst)
    return

def test_u(ceq):
    T = array([2500.0, 2200.0])
    Xs0 = array([[0.7, 0.2, 0.1],
                 [0.7, 0.2, 0.1]])
    print("Computing")
    u = ceq.batch_u(Xs0, T)
    print("Done: ", u)
    return

if __name__=='__main__':
    spnames = ['CO2', 'CO', 'O2']
    ceq = pyeq.EqCalculator(spnames)

    T = array([2500.0, 2000.0])
    p = array([0.1*101.35e3, 1.0*101.35e3])
    Xs0 = array([[1.0, 0.0, 0.0],
                 [1.0, 0.0, 0.0]])
    Xst = array([0.66108962603325838,0.22594024931116111,0.11297012465558055])

    test_pt(ceq, p, T, Xs0, Xst)
    print(' ')
    test_u(ceq)
