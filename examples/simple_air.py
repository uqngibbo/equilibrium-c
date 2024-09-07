"""
A simple air dissociation problem with eqc

@author: Nick Gibbons
"""
from numpy import array
import eqc

spnames = ['N2', 'O2', 'N', 'O', 'NO']
T = 2500.0
p = 0.1*101.35e3
Xs0 = array([0.767, 0.233, 0.0, 0.0, 0.0])

eq = eqc.EqCalculator(spnames)
Xs1 = eq.pt(p, T, Xs0, 1)

Xst = array([7.4785e-1, 2.0900e-1, 7.9320e-7, 2.0799e-2, 2.2349e-2])
Ys0 = eq.XtoY(Xs0)
print("Ys0: ", Ys0)
Ys1 = eq.XtoY(Xs1)
Yst = eq.XtoY(Xst)
print("Xs1: ", Xs1)

print("Mole Fractions:")
print("     " + ''.join(["      {:>2s}      ".format(spname) for spname in spnames]))
print("--------------------------------------------------------------------------")
print(" eqc:" + ''.join([" {:<10.8g}".format(Xs1s) for Xs1s in Xs1]))
print(" CEA:" + ''.join([" {:<10.8g}".format(Xs1s) for Xs1s in Xst]))

#print("Name  Mole Fraction   Target")
#for name, massf, target_massf in zip(spnames, Xs1, Xst):
#    print(" {:>2s}   {:<13.8g} - {:<13.8g}".format(name, massf, target_massf))

