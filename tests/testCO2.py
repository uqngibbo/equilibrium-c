"""
Automated test code for ceq

@author: Nick Gibbons
"""

from numpy import array, zeros, linspace
from pylab import semilogy, legend, xlabel, ylabel, title, show, ylim
from lewis_thermo import get_species
from os import system 
import pyeq
from re import search,split

def get_cea_test(p,T,spnames):
    with open('ceatests/template.inp') as fp: text=fp.read()
    text = text.replace('PRESSURE', str(p/101.35e3))
    text = text.replace('TEMPERATURE', str(T))
    with open('ceatests/co2.inp','w') as fp: fp.write(text)

    system('cea2 ceatests/co2') 
    with open('ceatests/co2.out') as fp: text=fp.read()

    regexstring = 'MOLE FRACTIONS\n\n( [\*]*[\S]+[\s]+[\s\.\d-]+\n)+\n'
    match = search(regexstring,text)
    if match==None: raise Exception("No match find in file!")

    lines = match.group().splitlines()
    keys = [i[:-8].strip().strip('*') for i in lines[2:-1]]
    values = [i[-8:].replace(' ','e').replace('-','e-') for i in lines[2:-1]]

    X = [0.0]*3
    for k,v in zip(keys,values):
        X[spnames.index(k)] = float(v)

    return X

def print_compare_table(X1, X2, name1, name2, T, spnames):
    print(' '*12 + ' '.join(['{:^27s}'.format(sp) for sp in spnames]))
    print(' '*12 + '  '.join(['{:^13s}{:^13s}'.format(name1, name2) for sp in spnames]))
    for Ti, X1i, X2i in zip(T, X1, X2):

        numbers = ' | '.join('{:.6e} {:.6e}'.format(*i) for i in zip(X1i,X2i))
        start = '{:7.2f} K : '.format(Ti)
        line = start + numbers
        print(line)

    print(" ")
    return


if __name__=='__main__':
    ps = [0.1*101.35e3, 1.0*101.35e3]
    Ts = linspace(301.0, 4000.0, 50)
    spnames = ['CO2', 'CO', 'O2']
    Xs0 = array([1.0, 0.0, 0.0])

    ceq = pyeq.EqCalculator(spnames)
    linecolours = ['red','blue','black']
    linestyles = ['-','--']
    markers = ['o','v']
    lines = []

    for i,p in enumerate(ps):
        print("p",p)
        pi = zeros(Ts.shape)
        pi[:] = p
        Xs0i = zeros((Ts.size, Xs0.size))
        Xs0i[:,:] = Xs0
        Xs1i = ceq.batch_pt(pi, Ts, Xs0i)

        icea = list(range(0,50,int(50/3)-1))
        Tcea = Ts[icea]
        Xcea = array([get_cea_test(p,Ti,spnames) for Ti in Tcea])
        Xceq = Xs1i[icea,:]
        print_compare_table(Xceq, Xcea, 'ceq','CEA', Tcea, spnames)

        plines = []
        pmarks = []
        for j in range(ceq.nsp):
            line = semilogy(Ts, Xs1i[:,j], linestyle=linestyles[i], color=linecolours[j])
            plines.extend(line)
            mark = semilogy(Tcea, Xcea[:,j], marker=markers[i],linestyle='None', color=linecolours[j])
            pmarks.extend(mark)
        lines.extend(plines+pmarks)

    lnames = ['']*ceq.nsp*3 + spnames
    ltitle = ' '+'           '.join(['{:8.2}'.format(p/101.35e3) for p in ps]) + '         (atm)\n'
    ltitle += ' ' + '    '.join(['ceq     CEA' for p in ps])
    xlabel('T (K)')
    ylabel('X')
    ylim(1e-7,10)
    legend(lines,lnames, loc='lower right', ncol=4, title=ltitle, columnspacing=0.0) 
    title('Carbon Dioxide Equilibrium Concentrations')
    show()
