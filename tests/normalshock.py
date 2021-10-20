"""
Application script: Compute a normal shock with ionization

@author: Nick Gibbons
"""

from numpy import zeros, array
from collections import namedtuple
from scipy.optimize import root, newton 
from copy import copy
import pyeq

Ru = 8.314
class GasState(object):
    def __init__(self, p, T, v, X0, ceq):
        p = max(100.0, p) # Optimizer sometimes likes to take WILD guesses on these
        T = max(20.0, T)  # which can go negative. This was enough to keep them converging.
        self.p=p; self.T=T; self.v=v; self.X0=X0; self.ceq=ceq
        
        X = ceq.pt(p, T, X0)
        h = ceq.get_h(X, T)
        s = ceq.get_s(X, T, p)
        cp= ceq.get_cp(X, T)
        Ms= ceq.M
        spnames = ceq.spnames

        Mmix = (Ms*X).sum()
        R = Ru/Mmix
        rho = p/R/T
        cv = cp-R
        k = cp/cv
        Y = X*Ms/Mmix

        self.X = X; self.h = h; self.s = s; self.cp=cp; self.Ms= Ms; self.spnames=spnames;
        self.Mmix = Mmix; self.R = R; self.rho = rho; self.cv=cv; self.k = k; self.Y = Y

        self.a = self.soundspeed()
        self.M = self.Mach_number()
        return

    def new_from_pTv(self, p, T, v):
        return GasState(p, T, v, X0=self.X0, ceq=self.ceq)

    def soundspeed(self):
        return (self.k*self.R*self.T)**0.5
    
    def Mach_number(self):
        a = self.soundspeed()
        return self.v/a

    def __repr__(self):
        s = [
        'p: {:12.3f} Pa  T: {:8.3f} K  rho: {:6.3f} kg/m3  v: {:8.3f} m/s'.format(self.p, self.T, self.rho, self.v),
        'h: {:12.3f} J/kg   s: {:8.3f} J/kg/K  Mmix: {:5.6f} kg/mol  M: {:5.6f}'.format(self.h, self.s, self.Mmix, self.M),
        ', '.join(['{}:{:8.7f}'.format(k,v) for k,v in zip(self.spnames, self.Y)])
        ]
        return '\n'.join(s)


def Rankine_Huginot_Error(s1, s2):
    mass = s1.rho*s1.v - s2.rho*s2.v
    mom = (s1.rho*s1.v**2 + s1.p)  - (s2.rho*s2.v**2 + s2.p)
    enth = (s1.h + s1.v**2/2.0)  - (s2.h + s2.v**2/2.0)
    error = array([mass**2, mom**2, enth**2])
    return error

def F(primitives, preshock):
    """ Compute the Rankine Huginot jump error from the preshock state in the incident shock frame """
    p,T,v = primitives
    postshock = preshock.new_from_pTv(p,T,v)
    error = Rankine_Huginot_Error(preshock, postshock) 
    return error 

def guess(s1):
    """ Guess initial postshock conditions using ideal gas behaviour """
    k = s1.k; R = s1.R; M1 = s1.M

    p2 = s1.p*(2*k*M1**2-(k-1))/(k+1)
    T2 = s1.T*(2*k*M1**2-(k-1))*((k-1)*M1**2+2.)/((k+1)**2*M1**2)
    M2= ((M1**2*(k-1)+2)/(2*k*M1**2-(k-1)))**0.5
    v2 = M2*(k*R*T2)**0.5
    return array([p2, T2, v2])

if __name__=='__main__':
    # T4 shot 12033 conditions
    spnames = ['N2', 'N2+', 'NO', 'NO+', 'O2', 'O2+', 'N', 'N+', 'O', 'O+', 'e-']
    X0 = array([0.77, 0.0, 0.0, 0.0, 0.23, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    X0/=(X0.sum())
    p1=40000; T1=300; vi=1777.73

    ceq = pyeq.EqCalculator(spnames)

    preshock = GasState(p=p1, T=T1, v=vi, X0=X0, ceq=ceq)

    # Compute a normal shock by solving for function F=0.0
    start = guess(preshock)
    sinfo = root(F, start, args=(preshock,))
    p,T,v = sinfo.x

    postshock = GasState(p=p, T=T, v=v, X0=X0, ceq=ceq)

    print("preshock:\n", preshock, '\n')
    print("postshock:\n", postshock, '\n')

    print("Mass fractions: ")
    for spname, mass_fraction in zip(postshock.spnames, postshock.Y):
        print(spname, mass_fraction)


