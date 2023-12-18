"""
Application script: Compute reflected shock tunnel output

@author: Nick Gibbons
"""

from numpy import zeros, array
from collections import namedtuple
from scipy.optimize import root, newton 
from copy import copy
from normalshock import *
import pyeq

def expand_to_a_ratio(pguess, s6, a_ratio):
    if pguess<1.0: p=10.0

    s7 = s6.expand_isentropically_to_p(pguess)
    rhoV_throat = s6.rho*s6.v
    rhoV_exit   = s7.rho*s7.v
    #print("called a ratio function with guess ", pguess, a_ratio)
    #print("Got result ", s7)
    #print("rhoV ", rhoV_exit, rhoV_throat)
    #print("")

    error = rhoV_exit/rhoV_throat - 1.0/a_ratio
    return error

def stna(p1, T1, vi, Ys1, pe, a_ratio):
    """ Shock tube problem with known nozzle area ratio """
    spnames = list(Ys1.keys())
    Y1 = array(list(Ys1.values()))

    ceq = pyeq.EqCalculator(spnames)

    Mmix = (Y1/ceq.M).sum()
    Mmix = Mmix**-1
    X0 = Y1*Mmix/ceq.M

    # Shock tube fill condition in incident shock frame (ISF)
    s1_isf = GasState.from_pTv(p=p1, T=T1, v=vi, X0=X0, ceq=ceq)

    # Compute state 2 in incident shock frame
    s2_isf = normal_shock(s1_isf)

    # Transform s2 into lab frame
    s2 = s2_isf.new_from_pTv(s2_isf.p, s2_isf.T, s1_isf.v - s2_isf.v)

    # Compute reflected shock state
    s5_rsf = reflected_normal_shock(s2)
    # Transform to lab frame
    s5 = copy(s5_rsf)
    s5.v = 0.0
    s5.M = 0.0

    # Compute relaxed stagnation state
    s5s_X, s5s_T = ceq.ps(pe, s5.s, X0)
    s5s = GasState.from_pTv(pe, s5s_T, 0.0, s5s_X, ceq=ceq)

    # Compute nozzle throat state
    expand_to_M1 = lambda p : (s5s.expand_isentropically_to_p(p)).M - 1.0
    p6 = newton(expand_to_M1, s5s.p/2.0)
    s6 = s5s.expand_isentropically_to_p(p6)

    # Compute nozzle exit state 
    a_ratio_error = lambda p : expand_to_a_ratio(p, s6, a_ratio)
    p7 = newton(a_ratio_error, s6.p*0.001)
    s7 = s6.expand_isentropically_to_p(p7)
    return [s1_isf, s2, s5, s5s, s6, s7]


def stnp(p1, T1, vi, Ys1, pe, pp_on_pe):
    spnames = list(Ys1.keys())
    Y1 = array(list(Ys1.values()))

    ceq = pyeq.EqCalculator(spnames)

    Mmix = (Y1/ceq.M).sum()
    Mmix = Mmix**-1
    X0 = Y1*Mmix/ceq.M

    # Shock tube fill condition in incident shock frame (ISF)
    s1_isf = GasState.from_pTv(p=p1, T=T1, v=vi, X0=X0, ceq=ceq)

    # Compute state 2 in incident shock frame
    s2_isf = normal_shock(s1_isf)

    # Transform s2 into lab frame
    s2 = s2_isf.new_from_pTv(s2_isf.p, s2_isf.T, s1_isf.v - s2_isf.v)

    # Compute reflected shock state
    s5_rsf = reflected_normal_shock(s2)
    # Transform to lab frame
    s5 = copy(s5_rsf)
    s5.v = 0.0
    s5.M = 0.0

    # Compute relaxed stagnation state
    s5s_X, s5s_T = ceq.ps(pe, s5.s, X0)
    s5s = GasState.from_pTv(pe, s5s_T, 0.0, s5s_X, ceq=ceq)

    # Compute nozzle throat state
    expand_to_M1 = lambda p : (s5s.expand_isentropically_to_p(p)).M - 1.0
    p6 = newton(expand_to_M1, s5s.p/2.0)
    s6 = s5s.expand_isentropically_to_p(p6)

    # Compute nozzle exit state 
    expand_to_pp_on_pe = lambda p : (s5s.expand_isentropically_to_p(p)).pitot_pressure()/s5s.p - pp_on_pe
    p7 = newton(expand_to_pp_on_pe, s5s.p*0.01)
    s7 = s5s.expand_isentropically_to_p(p7)
    return [s1_isf, s2, s5, s5s, s6, s7]

if __name__=='__main__':
    # T4 shot 12033 conditions
    #Ys0={'N2': 0.75518, 'Ar': 0.012916, 'CO2': 0.00048469, 'O2': 0.23142, 'NO':0.0, 'O':0.0}
    Ys1={'N2': 0.75518, 'Ar': 0.012916, 'CO2': 0.00048469, 'O2': 0.23142, 'NO':0.0}
    p1=40000; T1=300; vi=1777.73
    pe = 5.1687e+06
    pp_on_pe = 0.12
    states = stnp(p1, T1, vi, Ys1, pe, pp_on_pe)

    s1_isf, s2, s5, s5s, s6, s7 = states
    print("s1\n", s1_isf, '\n')
    print("s2\n", s2, '\n')
    print("s5\n", s5, '\n')
    print("s5s\n", s5s, '\n')
    print("s6\n", s6, '\n')
    print("s7\n", s7, '\n')

    #for sname,s in zip(['s1, s2, s5, s5s, s6, s7'], states):
    #    print(sname)
    #    print(s)
    #    print(" ")

