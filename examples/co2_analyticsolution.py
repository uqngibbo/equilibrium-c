"""
A simple CO2 dissociation problem, for the paper 

@author: Nick Gibbons
"""
from numpy import array,exp,cbrt,sqrt
import pyeq

#Gibbs free energy at standard pressure of one atmosphere
spnames =    [    "CO2",      "CO", "O2"]
gT = {1500 : [-396352.0, -243674.0, 0.0],
      2000 : [-396410.0, -285948.0, 0.0],
      2500 : [-396152.0, -327245.0, 0.0],
      3000 : [-395562.0, -367684.0, 0.0]}


nu = [-1.0,1.0, 0.5]
T = 2000.0
p = 0.1*101.35e3


Ru = 8.314
po = 101.35e3
G = sum(gj*nuj for gj,nuj in zip(gT[int(T)], nu))
Kp = exp(-G/Ru/T)

#alpha_coeffs = [(1.0 - p/po/Kp**2),
#                               0.0,
#                              -3.0,
#                               2.0]
#
#solutions = roots(alpha_coeffs)
#print("solutions", solutions)
#sol = solutions[-1]
#print("sol", sol)
#alpha = sol.real

A = (1.0 - p/po/Kp**2)
alpha = cbrt(sqrt(A**4-A**3) - A**2)/A + 1.0/cbrt(sqrt(A**4-A**3) - A**2)
print("sol: ", alpha)

total = 1.0 + alpha/2.0
Xans = [(1.0-alpha)/total, alpha/total, alpha/2.0/total]


ceq = pyeq.EqCalculator(spnames)
Xs0 = array([1.0, 0.0, 0.0])

Xs1 = ceq.pt(p, T, Xs0, 0)

print("Test: T={} K   p={} (Pa)".format(T, p))

print("Name  Mole Fraction")
for name, molef, molef2 in zip(spnames, Xs1, Xans):
    print(" {:>2s}   {:<13.8g}   {:<13.8g}".format(name,molef,molef2))

