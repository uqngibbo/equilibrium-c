"""
Application script: Compute shocktube problem data and compare to exp data. 

@author: Nick Gibbons
"""

import matplotlib.pyplot as plt
from numpy import pi
from shocktube import *
import eqc

# TODO Find the fix for the minus sign issues
plt.rcParams.update({'font.size': 12})
plt.rcParams['svg.fonttype'] = 'none'

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
r *= 1000.0

fig = plt.figure(figsize=(4.5,5))
axes = fig.gca()

rcfd = linspace(0.0, 100.0)

ppcfd_0 = 0.0*rcfd + pp_on_pe_0
ppcfd_6 = 0.0*rcfd + pp_on_pe_6

axes.errorbar(pp, r, xerr=pp*0.1, marker='.', color='black', linestyle='None', capsize=3, elinewidth=1, markeredgewidth=2,label = 'Shot 11310')
axes.plot(ppcfd_0, rcfd,'b-', linewidth=3, label="eqc (0mm)")
axes.plot(ppcfd_6, rcfd,'b--', linewidth=3, label="eqc (6mm)")
axes.fill_betweenx(rcfd, ppcfd_0, ppcfd_6, alpha=0.3, color='blue')

axes.set_ylabel('Radial Distance from Nozzle Axis (mm)')
axes.set_xlabel('Pressure \$p_7/p_5\$')
#axes.yaxis.set_label_position("right")
#axes.yaxis.tick_right()
#axes.set_title('Mach 7 Nozzle Validation')
axes.spines['top'].set_visible(False)
axes.spines['right'].set_visible(False)
axes.legend(framealpha=1.0, loc='lower left')
plt.tight_layout()
plt.grid()
plt.savefig('t4_11310_results.svg')
plt.show()
