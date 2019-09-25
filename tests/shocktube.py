"""
Application script: Compute reflected shock tunnel output

@author: Nick Gibbons
"""

from numpy import zeros, array
from collections import namedtuple
from scipy.optimize import root, newton 
from copy import copy
import pyeq

Ru = 8.314
class GasState(object):
    def __init__(self, p, T, v, X0):
        p = max(100.0, p) # Optimizer sometimes likes to take WILD guesses on these
        T = max(20.0, T)  # which can go negative. This was enough to keep them converging.
        self.p=p; self.T=T; self.v=v; self.X0=X0;
        
        X = self.get_X(p, T, X0)
        h = ceq.get_h(X, T)
        s = ceq.get_s(X, T, p)
        cp = ceq.get_cp(X, T)
        Mmix = (ceq.M*X).sum()

        R = Ru/Mmix
        rho = p/R/T
        cv = cp-R
        k = cp/cv
        Y = X*ceq.M/Mmix

        self.X = X; self.h = h; self.s = s; self.cp=cp; self.Mmix = Mmix;
        self.R = R; self.rho = rho; self.cv=cv; self.k = k; self.Y = Y

    def get_X(self, p, T, X0):
        return ceq.pt(p, T, X0)

    def __repr__(self):
        s = [
        'p: {:12.3f} Pa  T: {:8.3f} K  rho: {:6.3f} kg/m3  v: {:8.3f} m/s'.format(self.p, self.T, self.rho, self.v),
        'h: {:12.3f} J/kg   s: {:8.3f} J/kg/K  M: {:5.6f} kg/mol'.format(self.h, self.s, self.Mmix),
        ', '.join(['{}:{:8.7f}'.format(k,v) for k,v in zip(ceq.spnames, self.Y)])
        ]
        return '\n'.join(s)

class IdealGasState(GasState):
    def get_X(self, p, T, X0):
        return X0

def Rankine_Huginot_Error(s1, s2):
    mass = s1.rho*s1.v - s2.rho*s2.v
    mom = (s1.rho*s1.v**2 + s1.p)  - (s2.rho*s2.v**2 + s2.p)
    enth = (s1.h + s1.v**2/2.0)  - (s2.h + s2.v**2/2.0)
    error = array([mass**2, mom**2, enth**2])
    return error

def F(primitives, preshock):
    """ Compute the Rankine Huginot jump error from the preshock state in the incident shock frame """
    p,T,v = primitives
    postshock = GasState(p,T,v,preshock.X0)
    error = Rankine_Huginot_Error(preshock, postshock) 
    return error 

def F2(primitives, preshock):
    """ Compute the Rankine Huginot jump error from state 2 in the lab frame, assuming v5=0.0 """
    p,T,v = primitives # Note that this v is actually the shock velocity, not just v5_dash
    postshock_dash = GasState(p,T,v,preshock.X)
    
    preshock_dash = copy(preshock)
    preshock_dash.v += v
    error = Rankine_Huginot_Error(preshock_dash, postshock_dash) 
    return error

def guess(s1):
    """ Guess initial postshock conditions using ideal gas behaviour """
    k = s1.k 
    R = s1.R
    M1 = Mach_number(s1)

    p2 = s1.p*(2*k*M1**2-(k-1))/(k+1)
    T2 = s1.T*(2*k*M1**2-(k-1))*((k-1)*M1**2+2.)/((k+1)**2*M1**2)
    M2= ((M1**2*(k-1)+2)/(2*k*M1**2-(k-1)))**0.5
    v2 = M2*(k*R*T2)**0.5
    return array([p2, T2, v2])

def expand_isentropically_to_p(p, state):
    #state_s = ceq.get_s(state.X, state.T, state.p)
    state_ht = state.h + state.v**2/2.0
    X, T = ceq.ps(p, state.s, state.X)
    h = ceq.get_h(X, T)
    if state_ht<h: raise Exception("Thermo Error: Too much expansion requests: p={}".format(p))

    v = (2.0*(state_ht - h))**0.5
    newstate = GasState(p, T, v, X)
    return newstate

def soundspeed(state):
    return (state.k*state.R*state.T)**0.5

def Mach_number(state):
    a = soundspeed(state)
    return state.v/a

def pitot_pressure(state):
    M = Mach_number(state)
    if M>1.0: # FIXME: Replace with a single normal shock function
        start = guess(state)
        sinfo = root(F, start, args=(state,))
        p,T,v = sinfo.x
        postshock = GasState(p=p, T=T, v=v, X0=state.X)
    else:
        postshock = copy(state)

    ht = postshock.h + postshock.v**2/2.0
    sps = ceq.get_s(postshock.X, postshock.T, postshock.p)
    stagnation_enthalpy_error = lambda p : ceq.get_h(*ceq.ps(p, sps, postshock.X)) - ht
    pstag = newton(stagnation_enthalpy_error, postshock.p*1.1) 
    return pstag


#Ys0={'N2': 0.75518, 'Ar': 0.012916, 'CO2': 0.00048469, 'O2': 0.23142, 'NO':0.0, 'O':0.0}
Ys0={'N2': 0.75518, 'Ar': 0.012916, 'CO2': 0.00048469, 'O2': 0.23142, 'NO':0.0}
spnames = list(Ys0.keys())
Y0 = array(list(Ys0.values()))
print(spnames)

ceq = pyeq.EqCalculator(spnames)

Mmix = (Y0/ceq.M).sum()
Mmix = Mmix**-1
X0 = Y0*Mmix/ceq.M
print(X0)

# Shock tube fill condition in incident shock frame (ISF)
s1_isf = GasState(p=40000, T=300, v=1777.73, X0=X0)
print("State 1: Pre shock condition")
print(s1_isf)
print(" ")


# Compute state 2 in incident shock frame
start = guess(s1_isf)
print("Incident shock initial guess", start)

sinfo = root(F, start, args=(s1_isf,))
p,T,v = sinfo.x
print("Finish", p,T,v)

s2_isf = GasState(p=p, T=T, v=v, X0=X0)
print("s2_isf")
print(s2_isf)
print(" ")

# Transform s2 into lab frame
s2 = copy(s2_isf)
s2.v = s1_isf.v - s2_isf.v
print("s2")
print(s2)
print(" ")

# Compute reflected shock state
print("")
print("Computing reflected shock")
start2 = guess(s2) # Sketchy guess that ignores the reflected velocity 
sinfo = root(F2, start2, args=(s2,))
p,T,v = sinfo.x
s5_dash = GasState(p=p, T=T, v=v, X0=X0)

s5 = copy(s5_dash)
s5.v = 0.0
print("State 5: Post reflection shock condition")
print(s5)
print("Shock Velocity: ", s5_dash.v)
print(" ")

# Compute relaxed stagnation state
pe = 5.1687e+06
s5_entropy = ceq.get_s(s5.X, s5.T, s5.p)
s5s_X, s5s_T = ceq.ps(pe, s5_entropy, X0)
s5s = GasState(pe, s5s_T, 0.0, s5s_X)
print("State 5s: Relax to pe")
print(s5s)
print(" ")

# Compute nozzle throat state
expand_to_M1 = lambda p : Mach_number(expand_isentropically_to_p(p, s5s)) - 1.0
p6 = newton(expand_to_M1, s5s.p/2.0)
s6 = expand_isentropically_to_p(p6, s5s)
print("State 6: Nozzle throat")
print(s6)
print("M:", Mach_number(s6))
print(" ")

# Compute nozzle exit state # FIXME: non cheat version needs to solve for pp_on_pe
pp_on_pe = 0.12
expand_to_pp_on_pe = lambda p : pitot_pressure(expand_isentropically_to_p(p, s5s))/s5s.p - pp_on_pe
p7 = newton(expand_to_pp_on_pe, s5s.p*0.01)
#p7 = 32841.0
s7 = expand_isentropically_to_p(p7, s5s)
M7 = Mach_number(s7)
pp7 = pitot_pressure(s7)
print("State 7: Nozzle Output")
print(s7)
print("M:", M7)
print("pp7", pp7, pp7/s5s.p)
print(" ")

