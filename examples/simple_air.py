"""
A simple air dissociation problem with ceq

@author: Nick Gibbons
"""
from numpy import array
import pyeq

spnames = ['N2', 'O2', 'N', 'O', 'NO']
T = 2500.0
p = 0.1*101.35e3
Xs0 = array([0.767, 0.233, 0.0, 0.0, 0.0])

ceq = pyeq.EqCalculator(spnames)
Xs1 = ceq.pt(p, T, Xs0, 0)

Xst = array([7.4785e-1, 2.0900e-1, 7.9320e-7, 2.0799e-2, 2.2349e-2])
Ys0 = ceq.XtoY(Xs0)
print("Ys0: ", Ys0)
Ys1 = ceq.XtoY(Xs1)
Yst = ceq.XtoY(Xst)
print("Xs1: ", Xs1)

print("Mole Fractions:")
print("     " + ''.join(["      {:>2s}      ".format(spname) for spname in spnames]))
print("--------------------------------------------------------------------------")
print(" ceq:" + ''.join([" {:<10.8g}".format(Xs1s) for Xs1s in Xs1]))
print(" CEA:" + ''.join([" {:<10.8g}".format(Xs1s) for Xs1s in Xst]))

#print("Name  Mole Fraction   Target")
#for name, massf, target_massf in zip(spnames, Xs1, Xst):
#    print(" {:>2s}   {:<13.8g} - {:<13.8g}".format(name, massf, target_massf))

