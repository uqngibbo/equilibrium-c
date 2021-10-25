"""
Application script: Compute a normal shock with ionization

v2: With pressure error method.

@author: Nick Gibbons
"""

from numpy import zeros, array, linspace
from collections import namedtuple
from scipy.optimize import newton 
import pyeq

Ru = 8.314
class GasState(object):
    def __init__(self, rho, T, p, u, v, X, ceq):
        h = ceq.get_h(X, T)
        s = ceq.get_s(X, T, p)
        cp= ceq.get_cp(X, T)
        Ms= ceq.M
        spnames = ceq.spnames

        Mmix = (Ms*X).sum()
        R = Ru/Mmix
        cv = cp-R
        k = cp/cv
        Y = X*Ms/Mmix

        self.T = T; self.p = p; self.u = u; self.v = v; self.X = X; self.ceq = ceq;
        self.h = h; self.s = s; self.cp=cp; self.Ms= Ms; self.spnames=spnames;
        self.Mmix = Mmix; self.R = R; self.rho = rho; self.cv=cv; self.k = k; self.Y = Y; self.X0=X;

        self.a = self.soundspeed()
        self.M = self.Mach_number()

    @classmethod
    def from_pTv(cls, p, T, v, X0, ceq):
        p = max(1.0, p)  # Optimizer sometimes likes to take WILD guesses on these
        T = max(1.0, T)  # which can go negative. This was enough to keep them converging.
        
        X = ceq.pt(p, T, X0)
        u = ceq.get_u(X, T)

        Mmix = (ceq.M*X).sum()
        R = Ru/Mmix
        rho = p/R/T

        return cls(rho, T, p, u, v, X, ceq)

    @classmethod
    def from_rhouv(cls, rho, u, v, X0, ceq):
        X, T = preshock.ceq.rhou(rho, u, X0, 0)
        Mmix = (ceq.M*X).sum()
        R = Ru/Mmix
        p = rho*R*T
        return cls(rho, T, p, u, v, X, ceq)

    def new_from_pTv(self, p, T, v):
        return GasState.from_pTv(p, T, v, X0=self.X0, ceq=self.ceq)

    def new_from_rhouv(self, rho, u, v):
        return GasState.from_rhouv(rho, u, v, X0=self.X0, ceq=self.ceq)

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

def fluxes(gs):
    rho = gs.rho; v = gs.v; p = gs.p;

    Fmass = rho*v
    Fspecies = [Fmass*Y for Y in gs.Y]
    Fmom = rho*v**2 + p
    Fenergy = v*(rho*gs.u + 0.5*rho*v**2 + p)
    return Fmass, Fspecies, Fmom, Fenergy

def rhopu_from_v(v, preshock):
    Fmass, Fspecies, Fmom, Fenergy = fluxes(preshock)
    Y0 = array([Fsps/Fmass for Fsps in Fspecies])
    X0 = preshock.ceq.YtoX(Y0)

    rho = Fmass/v
    p = Fmom - rho*v**2
    u = (Fenergy - 0.5*rho*v**3 - p*v)/(Fmass)
    return rho, p, u

def pressure_error_from_v(v, preshock):
    v = max(10.0, v) # We have to cap this at a positive number, to prevent negative rho
    rho, p, u = rhopu_from_v(v, preshock)
    postshock = preshock.new_from_rhouv(rho, u, v)

    # Here is where things start to get weird...
    
    # We've computed the entire gas state using the fluxes and a guessed velocity
    # v. But that makes no sense because *any* v that we chose would have given us
    # an answer??? The key is that there is one equation we haven't used yet, the
    # ideal gas equation of state p = rho*R*T. Inside that new_from_rhouv call
    # above, we computed a pressure using the ideal gas equation.  That pressure is
    # actually inconsistent with the pressure calculated above from the momentum flux.
    # And so: We can adjust the guess velocity v until the two match up and that
    # will give the correct postshock state.

    return (p - postshock.p)/p

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
    p1=40000; T1=300; vi=3000.0

    ceq = pyeq.EqCalculator(spnames)

    preshock = GasState.from_pTv(p=p1, T=T1, v=vi, X0=X0, ceq=ceq)
    start = guess(preshock)
    print("preshock:\n", preshock, '\n')
    print("guess:\n p={} T={} v={}\n".format(*start))

    # Compute a normal shock by solving for function f=0.0
    f = lambda v : pressure_error_from_v(v, preshock)
    vpostshock = newton(f, start[2]) # guess velocity is start[2]

    # With the correct velocity found, compute the postshock state
    rho, p, u = rhopu_from_v(vpostshock, preshock)
    postshock = preshock.new_from_rhouv(rho, u, vpostshock)
    ionisation_fraction = postshock.X[spnames.index('e-')]*2

    print("postshock:\n", postshock, '\n')

    print("Ionization fraction: {:g} %".format(ionisation_fraction*100))



