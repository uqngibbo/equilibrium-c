"""
Application script: Compute shocktube problem data and compare to exp data. 

@author: Nick Gibbons
"""

from numpy import pi
from shocktube import *
import pyeq


# T4 shot 11310 conditions, from M7_pitot_survey, table 1
Ys1={'N2': 0.767, 'O2': 0.233, 'NO':0.0, 'N':0.0, 'O':0.0, 'NO':0.0}
p1=200e3; T1=300; vi=1669.0
pe = 19.53e6
pp_on_pe = 0.107
r_throat = 1.050000e-2 # From t4m7 eilmer example
r_exit   = 1.365830e-1 - 0e-3 # From t4m7 eilmer example, minus wilson's estimated 6mm displacement thickness
a_ratio = (pi*r_exit**2)/(pi*r_throat**2)
#states = stnp(p1, T1, vi, Ys1, pe, pp_on_pe)
states = stna(p1, T1, vi, Ys1, pe, a_ratio)

s1_isf, s2, s5, s5s, s6, s7 = states
print("s1\n", s1_isf, '\n')
print("s2\n", s2, '\n')
print("s5\n", s5, '\n')
print("s5s\n", s5s, '\n')
print("s6\n", s6, '\n')
print("s7\n", s7, '\n')

pitot = s7.pitot_pressure()
print("pitot pressure: ", pitot)

print("pitot/stag ratio: ", pitot/pe)
