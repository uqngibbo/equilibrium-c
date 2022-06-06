"""
A simple air dissociation problem with ceq

@author: Nick Gibbons
"""
from numpy import array
import pyeq

spnames = ['N2', 'O2', 'N', 'O', 'NO']
T = 2500.0
p = 0.1*101.35e3
Xs0 = array([0.76, 0.23, 0.0, 0.0, 0.0])
Xst = array([7.48543073e-01, 2.08366049e-01, 7.93468988e-07, 2.07645979e-02, 2.23254866e-02])

ceq = pyeq.EqCalculator(spnames)

Xs1 = ceq.pt(p, T, Xs0, 0)
Ys0 = ceq.XtoY(Xs0)
Ys1 = ceq.XtoY(Xs1)
Yst = ceq.XtoY(Xst)

print("Name  Mass Fraction   Target")
for name, massf, target_massf in zip(spnames, Ys1, Ys1):
    print(" {:>2s}   {:<13.8g} - {:<13.8g}".format(name, massf, target_massf))

