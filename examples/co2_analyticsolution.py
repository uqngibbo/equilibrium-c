"""
A simple CO2 dissociation problem, from:
"An Introduction to Combustion"
Stephen R. Turns, Third Edition
Example 2.7

@author: Nick Gibbons
"""
from numpy import array,exp,cbrt,roots,linspace,zeros,log10
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import eqc

# TODO Find the fix for the minus sign issues
plt.rcParams.update({'font.size': 12})
plt.rcParams['svg.fonttype'] = 'none'

def make_2D_legend(axes, col_labels, row_labels, lines):
    row_names = ['']*(len(col_labels)-1)*len(row_labels) + row_labels
    spacing = max(len(label) for label in col_labels)
    col_names = '     '.join(col_labels)

    print(row_names)
    legend = axes.legend(lines, row_names, loc='lower left', ncol=len(col_labels), title=col_names, columnspacing=1.0, framealpha=1.0, handlelength=0.49*spacing)
    legend._legend_box.align = "left"
    return legend

sqrt = np.emath.sqrt

# Gibbs energy at standard pressure of one atmosphere, from Turns' Appendix A
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

lines = []
markers=[]
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

        # Python has trouble with the complex numbers here. Let's just use the 
        # roots function from numpy 
        #alpha = cbrt(sqrt(A**4-A**3) - A**2)/A + 1.0/cbrt(sqrt(A**4-A**3) - A**2)
        solutions = roots(alpha_coeffs)
        alpha = solutions[-1].real

        total = 1.0 + alpha/2.0
        Xans = [(1.0-alpha)/total, alpha/total, alpha/2.0/total]
        Xss.append(Xans)

        print("Test: T={} K   p={} (Pa)".format(T, p))


    Xss = array(Xss)
    markers.extend(axes[0].semilogy(Ts, Xss[:,0], linestyle="none", marker='o', color='black', markersize=4.0+1.0*i))
    axes[1].semilogy(Ts, Xss[:,1], linestyle="none", marker='o', color='blue', markersize=4.0+1.0*i)
    axes[2].semilogy(Ts, Xss[:,2], linestyle="none", marker='o', color='red', markersize=4.0+1.0*i)

    Tss = linspace(Ts[0], Ts[-1])
    pss = Tss*0.0 + p
    Xs0 = zeros((Tss.size, len(spnames)))
    Xs0[:,0] = 1.0
    eq = eqc.EqCalculator(spnames)
    Xs1 = eq.batch_pt(pss, Tss, Xs0, 0)
    lines.extend(axes[0].semilogy(Tss, Xs1[:,0], color='black', linewidth=0.5+0.5*i, label="{:3.1f} Atms".format(p/po)))
    axes[1].semilogy(Tss, Xs1[:,1], color='blue', linewidth=0.5+0.5*i)
    axes[2].semilogy(Tss, Xs1[:,2], color='red', linewidth=0.5+0.5*i)

symbols = markers+lines
make_2D_legend(axes[0], ['Turns', 'eq-c'], ['0.1 Atm', '1 Atm', '10 Atm','100 Atm'], symbols)

axes[0].set_title('CO2')
axes[1].set_title('CO')
axes[2].set_title('O2')

axes[0].set_ylabel('Mole Fraction')
axes[0].set_xlabel('Temperature (K)')
axes[1].set_xlabel('Temperature (K)')
axes[2].set_xlabel('Temperature (K)')

axes[0].set_ylim(0.09, 1.1)
axes[1].set_ylim(1e-5, 1.0)
axes[2].set_ylim(1e-5, 1.0)

fig.tight_layout()

savefig = True
if savefig:
    for ax in axes:
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '\$10^{{{:d}}}\$'.format(int(log10(y)))))
        for tick in ax.yaxis.get_major_ticks(): tick.label.set_fontsize(4)
    plt.savefig('co2.svg')

plt.show()
