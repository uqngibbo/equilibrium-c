"""
Application script: Compute reflected shock tunnel output

@author: Nick Gibbons
"""

from numpy import zeros, array
from collections import namedtuple
from scipy.optimize import root, newton 
from copy import copy
import eqc

Ru = 8.314
class GasState(object):
    def __init__(self, p, T, v, X0, eq):
        p = max(100.0, p) # Optimizer sometimes likes to take WILD guesses on these
        T = max(20.0, T)  # which can go negative. This was enough to keep them converging.
        self.p=p; self.T=T; self.v=v; self.X0=X0; self.eq=eq
        
        X = eq.pt(p, T, X0)
        h = eq.get_h(X, T)
        s = eq.get_s(X, T, p)
        cp= eq.get_cp(X, T)
        Ms= eq.M
        spnames = eq.spnames

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
        return GasState(p, T, v, X0=self.X0, eq=self.eq)

    def soundspeed(self):
        return (self.k*self.R*self.T)**0.5
    
    def Mach_number(self):
        a = self.soundspeed()
        return self.v/a

    def expand_isentropically_to_p(self, p):
        state_ht = self.h + self.v**2/2.0
        X, T = self.eq.ps(p, self.s, self.X0)
        h = self.eq.get_h(X, T)
        if state_ht<h: raise Exception("Thermo Error: Too much expansion requested: p={}".format(p))

        v = (2.0*(state_ht - h))**0.5
        newstate = self.new_from_pTv(p, T, v)
        return newstate
    
    def pitot_pressure(self):
        if self.M>1.0: # FIXME: Replace with a single normal shock function
            start = guess(self)
            sinfo = root(F, start, args=(self,))
            p,T,v = sinfo.x
            postshock = self.new_from_pTv(p=p, T=T, v=v)
        else:
            postshock = copy(self)
    
        ht = postshock.h + postshock.v**2/2.0
        sps = postshock.s
        stagnation_enthalpy_error = lambda p : self.eq.get_h(*self.eq.ps(p, sps, postshock.X0)) - ht
        pstag = newton(stagnation_enthalpy_error, postshock.p*1.1) 
        return pstag

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

def F2(primitives, preshock):
    """ Compute the Rankine Huginot jump error from state 2 in the lab frame, assuming v5=0.0 """
    p,T,v = primitives # Note that this v is actually the shock velocity, not just v5_dash
    postshock_dash = preshock.new_from_pTv(p,T,v)
    
    preshock_dash = copy(preshock)
    preshock_dash.v += v
    error = Rankine_Huginot_Error(preshock_dash, postshock_dash) 
    return error

def guess(s1):
    """ Guess initial postshock conditions using ideal gas behaviour """
    k = s1.k; R = s1.R; M1 = s1.M

    p2 = s1.p*(2*k*M1**2-(k-1))/(k+1)
    T2 = s1.T*(2*k*M1**2-(k-1))*((k-1)*M1**2+2.)/((k+1)**2*M1**2)
    M2= ((M1**2*(k-1)+2)/(2*k*M1**2-(k-1)))**0.5
    v2 = M2*(k*R*T2)**0.5
    return array([p2, T2, v2])

def stnp(p1, T1, vi, Ys1, pe, pp_on_pe):
    spnames = list(Ys1.keys())
    Y1 = array(list(Ys1.values()))
    print(spnames)

    eq = eqc.EqCalculator(spnames)

    Mmix = (Y1/eq.M).sum()
    Mmix = Mmix**-1
    X0 = Y1*Mmix/eq.M
    print(X0)

    # Shock tube fill condition in incident shock frame (ISF)
    s1_isf = GasState(p=p1, T=T1, v=vi, X0=X0, eq=eq)
    print("s1\n", s1_isf, '\n')

    # Compute state 2 in incident shock frame
    start = guess(s1_isf)
    sinfo = root(F, start, args=(s1_isf,))
    p,T,v = sinfo.x
    s2_isf = GasState(p=p, T=T, v=v, X0=X0, eq=eq)
    print("s2_isf\n", s2_isf, '\n')

    # Transform s2 into lab frame
    s2 = copy(s2_isf)
    s2 = s2_isf.new_from_pTv(s2.p, s2.T, s1_isf.v - s2_isf.v)
    print("s2\n", s2, '\n')

    # Compute reflected shock state
    print("s2.M", s2.M)
    start2 = guess(s2) # Sketchy guess that ignores the reflected velocity 
    print('start2 guess', start2)
    sinfo = root(F2, start2, args=(s2,))
    p,T,v = sinfo.x
    print("F2 pTv", p,T,v)
    s5_dash = GasState(p=p, T=T, v=v, X0=X0, eq=eq)

    s5 = copy(s5_dash)
    s5.v = 0.0
    s5.M = 0.0
    print("s5\n", s5, '\n')

    # Compute relaxed stagnation state
    s5s_X, s5s_T = eq.ps(pe, s5.s, X0)
    s5s = GasState(pe, s5s_T, 0.0, s5s_X, eq=eq)
    print("s5s\n", s5s, '\n')

    # Compute nozzle throat state
    expand_to_M1 = lambda p : (s5s.expand_isentropically_to_p(p)).M - 1.0
    p6 = newton(expand_to_M1, s5s.p/2.0)
    s6 = s5s.expand_isentropically_to_p(p6)
    print("s6\n", s6, '\n')

    # Compute nozzle exit state 
    expand_to_pp_on_pe = lambda p : (s5s.expand_isentropically_to_p(p)).pitot_pressure()/s5s.p - pp_on_pe
    p7 = newton(expand_to_pp_on_pe, s5s.p*0.01)
    s7 = s5s.expand_isentropically_to_p(p7)
    print("s7\n", s7, '\n')
    return [s1_isf, s2, s5, s5s, s6, s7]

if __name__=='__main__':
    # T4 shot 12033 conditions
    #Ys0={'N2': 0.75518, 'Ar': 0.012916, 'CO2': 0.00048469, 'O2': 0.23142, 'NO':0.0, 'O':0.0}
    Ys1={'N2': 0.75518, 'Ar': 0.012916, 'CO2': 0.00048469, 'O2': 0.23142, 'NO':0.0}
    p1=40000; T1=300; vi=1777.73
    pe = 5.1687e+06
    pp_on_pe = 0.12
    states = stnp(p1, T1, vi, Ys1, pe, pp_on_pe)

    #for sname,s in zip(['s1, s2, s5, s5s, s6, s7'], states):
    #    print(sname)
    #    print(s)
    #    print(" ")

