"""
Application script: Compute shocktube problem data and compare to exp data. 

@author: Nick Gibbons
"""

import matplotlib.pyplot as plt
from numpy import pi
from shocktube import *
import pyeq


# T4 shot 11310 conditions, from M7_pitot_survey, table 1
Ys1={'N2': 0.767, 'O2': 0.233, 'NO':0.0, 'N':0.0, 'O':0.0, 'NO':0.0}
p1=200e3; T1=300; vi=1669.0
pe = 19.53e6
r_throat = 1.050000e-2 # From t4m7 eilmer example
r_exit   = 1.365830e-1 - 6e-3 # From t4m7 eilmer example, minus wilson's estimated 6mm displacement thickness
a_ratio = (pi*r_exit**2)/(pi*r_throat**2)
states = stna(p1, T1, vi, Ys1, pe, a_ratio)

s1_isf, s2, s5, s5s, s6, s7 = states

pitot_6 = s7.pitot_pressure()
pp_on_pe_6 = pitot_6/pe
print("pitot pressure: ", pitot_6)
print("pitot/stag ratio: ", pp_on_pe_6)

r_throat = 1.050000e-2 # From t4m7 eilmer example
r_exit   = 1.365830e-1 - 0e-3 # From t4m7 eilmer example, with no correction for boundary layer
a_ratio = (pi*r_exit**2)/(pi*r_throat**2)
states = stna(p1, T1, vi, Ys1, pe, a_ratio)

s1_isf, s2, s5, s5s, s6, s7 = states

pitot_0 = s7.pitot_pressure()
pp_on_pe_0 = pitot_0/pe
print("pitot pressure: ", pitot_0)
print("pitot/stag ratio: ", pp_on_pe_0)

with open('t4_11310_pitot_data.csv') as fp:
    lines = [line.strip().split(',') for line in fp]

numbers = [list(map(float,line)) for line in lines]
r, pp = [array(i) for i in zip(*numbers)]

fig = plt.figure()
axes = fig.gca()

rcfd = linspace(0.0, 0.1)
ppcfd_0 = 0.0*rcfd + pp_on_pe_0
ppcfd_6 = 0.0*rcfd + pp_on_pe_6
axes.errorbar(r, pp, yerr=pp*0.1, fmt='k.', capsize=3, elinewidth=1, markeredgewidth=2,label = 'Shot 11310')
axes.plot(rcfd, ppcfd_0, 'b-', linewidth=2, label="ceq (0mm)")
axes.plot(rcfd, ppcfd_6, 'r-', linewidth=2, label="ceq (6mm)")
axes.set_xlabel('Radial Distance from Nozzle Axis (m)')
axes.set_ylabel('pitot pressure/stagnation pressure')
axes.set_title('Mach 7 Nozzle Validation')
axes.legend(framealpha=1.0)
plt.grid()
plt.show()
