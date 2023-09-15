"""
A simple CO2 dissociation problem, for the paper 

@author: Nick Gibbons
"""
from numpy import array,exp,cbrt,roots,linspace,zeros
import numpy as np
import matplotlib.pyplot as plt
import pyeq

sqrt = np.emath.sqrt
#Gibbs free energy at standard pressure of one atmosphere
spnames =    [    "CO2",      "CO", "O2"]
gT = {1500 : [-396352.0, -243674.0, 0.0],
      2000 : [-396410.0, -285948.0, 0.0],
      2500 : [-396152.0, -327245.0, 0.0],
      3000 : [-395562.0, -367684.0, 0.0]}

Ru = 8.314
po = 101.35e3
nu = [-1.0,1.0, 0.5]
ps= [0.1*101.35e3, 1*101.35e3, 10*101.35e3, 100*101.35e3]

fig = plt.figure(figsize=(12,4))
axes = fig.subplots(1,3)

for i,p in enumerate(ps):
    Ts = sorted(gT.keys())
    Xss = []
    for T in Ts:
        G = sum(gj*nuj for gj,nuj in zip(gT[T], nu))
        Kp = exp(-G/Ru/T)

        alpha_coeffs = [
            (1.0 - p/po/Kp**2),
            0.0,
           -3.0,
            2.0]

        solutions = roots(alpha_coeffs)
        alpha = solutions[-1].real
        #alpha = cbrt(sqrt(A**4-A**3) - A**2)/A + 1.0/cbrt(sqrt(A**4-A**3) - A**2)

        total = 1.0 + alpha/2.0
        Xans = [(1.0-alpha)/total, alpha/total, alpha/2.0/total]
        Xss.append(Xans)

        print("Test: T={} K   p={} (Pa)".format(T, p))


    Xss = array(Xss)
    axes[0].semilogy(Ts, Xss[:,0], linestyle="none", marker='o', color='black', linewidth=0.5+0.5*i)
    axes[1].semilogy(Ts, Xss[:,1], linestyle="none", marker='o', color='blue', linewidth=0.5+0.5*i)
    axes[2].semilogy(Ts, Xss[:,2], linestyle="none", marker='o', color='red', linewidth=0.5+0.5*i)

    Tss = linspace(Ts[0], Ts[-1])
    pss = Tss*0.0 + p
    Xs0 = zeros((Tss.size, len(spnames)))
    Xs0[:,0] = 1.0
    ceq = pyeq.EqCalculator(spnames)
    Xs1 = ceq.batch_pt(pss, Tss, Xs0, 0)
    axes[0].semilogy(Tss, Xs1[:,0], color='black', linewidth=0.5+0.5*i, label="{:3.1f} Atms".format(p/po))
    axes[1].semilogy(Tss, Xs1[:,1], color='blue', linewidth=0.5+0.5*i)
    axes[2].semilogy(Tss, Xs1[:,2], color='red', linewidth=0.5+0.5*i)

axes[0].set_title('CO2')
axes[1].set_title('CO')
axes[2].set_title('O2')

axes[0].legend(loc=0)
axes[0].set_ylabel('Mole Fraction')
axes[0].set_xlabel('Temperature (K)')
axes[1].set_xlabel('Temperature (K)')
axes[2].set_xlabel('Temperature (K)')
fig.tight_layout()
plt.show()
