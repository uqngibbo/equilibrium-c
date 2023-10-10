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
        X, T = ceq.rhou(rho, u, X0, 0)
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

    def expand_isentropically_to_p(self, p):
        state_ht = self.h + self.v**2/2.0
        X, T = self.ceq.ps(p, self.s, self.X0)
        h = self.ceq.get_h(X, T)
        if state_ht<h: raise Exception("Thermo Error: Too much expansion requested: p={}".format(p))

        v = (2.0*(state_ht - h))**0.5
        newstate = self.new_from_pTv(p, T, v)
        return newstate
    
    def pitot_pressure(self):
        if self.M>1.0:
            postshock = normal_shock(self)
        else:
            postshock = copy(self)
    
        ht = postshock.h + postshock.v**2/2.0
        sps = postshock.s
        stagnation_enthalpy_error = lambda p : self.ceq.get_h(*self.ceq.ps(p, sps, postshock.X0)) - ht
        pstag = newton(stagnation_enthalpy_error, postshock.p*1.1) 
        return pstag

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
    #Y0 = array([Fsps/Fmass for Fsps in Fspecies])
    #X0 = preshock.ceq.YtoX(Y0)

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

def pressure_error_from_v_reflected(v_reflected, preshock):
    """ For a reflected shock, we know that the postshock gas should be stationary """
    # This v is the velocity of the incoming flow in the reflected shock reference frame
    v2_dash = preshock.v + v_reflected
    v2_dash = max(1.0, v2_dash)
    preshock_rsf = preshock.new_from_pTv(preshock.p, preshock.T, v2_dash)

    # In the reflected shock frame, the postshock gas is moving at v_reflected
    rho, p, u = rhopu_from_v(v_reflected, preshock_rsf)
    postshock = preshock_rsf.new_from_rhouv(rho, u, v_reflected)

    # As above, the difference between the poshock pressure and p is an indication
    # of how good our guess of v_reflected is
    return (p - postshock.p)/p

def guess(s1):
    """ Guess initial postshock conditions using ideal gas behaviour """
    k = s1.k; R = s1.R; M1 = s1.M

    p2 = s1.p*(2*k*M1**2-(k-1))/(k+1)
    T2 = s1.T*(2*k*M1**2-(k-1))*((k-1)*M1**2+2.)/((k+1)**2*M1**2)
    M2= ((M1**2*(k-1)+2)/(2*k*M1**2-(k-1)))**0.5
    v2 = M2*(k*R*T2)**0.5
    return array([p2, T2, v2])

def guess_reflected(s1):
    """ Guess post reflected shock conditions using ideal gas behaviour """
    k = s1.k; R = s1.R; M2 = s1.M;

    # I got this expression from Ingo Jahn's Drummon Tunnel prac sheet
    # It's the velocity of the incoming gas in the reflected shock 
    # reference frame, assuming M2 is the post incident shock velocity 
    # in the Lab frame. Okay.
    Mr = (k+1)/4.0*M2 + (1+(k+1)**2/16.0*M2**2)**0.5

    p5 = s1.p*(2*k*Mr**2-(k-1))/(k+1)
    T5 = s1.T*(2*k*Mr**2-(k-1))*((k-1)*Mr**2+2.)/((k+1)**2*Mr**2)
    M5_dash = ((Mr**2*(k-1)+2)/(2*k*Mr**2-(k-1)))**0.5
    v5_dash = M5_dash*(k*R*T5)**0.5
    return array([p5, T5, v5_dash])

def normal_shock(preshock, verbose=False):
    start = guess(preshock)
    if verbose:
        print("guess:\n p={} T={} v={}\n".format(*start))

    # Compute a normal shock by solving for function f=0.0
    f = lambda v : pressure_error_from_v(v, preshock)
    vpostshock = newton(f, start[2]/2) # guess velocity is start[2]

    # With the correct velocity found, compute the postshock state
    rho, p, u = rhopu_from_v(vpostshock, preshock)
    postshock = preshock.new_from_rhouv(rho, u, vpostshock)
    return postshock

def reflected_normal_shock(preshock, verbose=False):
    """
    Do a normal shock reflection where the post shock gas
    is stagnant in the lab reference frame.
    """
    start = guess_reflected(preshock)
    if verbose:
        print("guess:\n p={} T={} v={}\n".format(*start))

    # Compute a normal shock by solving for function f=0.0
    f = lambda v_reflected : pressure_error_from_v_reflected(v_reflected, preshock)
    v_reflected = newton(f, start[2]) # guess velocity is start[2]

    # With the correct velocity found, compute the postshock state
    v2_dash = preshock.v + v_reflected
    preshock_rsf = preshock.new_from_pTv(preshock.p, preshock.T, v2_dash)
    rho, p, u = rhopu_from_v(v_reflected, preshock_rsf)
    postshock = preshock_rsf.new_from_rhouv(rho, u, v_reflected)

    return postshock


if __name__=='__main__':
    # Arbitrary conditions
    spnames = ['N2', 'N2+', 'NO', 'NO+', 'O2', 'O-', 'N', 'N+', 'O', 'O+', 'e-']
    X0 = array([0.767, 0.0, 0.0, 0.0, 0.233, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    X0/=(X0.sum())
    p1=600; T1=300; vi=6000.0

    ceq = pyeq.EqCalculator(spnames)

    preshock = GasState.from_pTv(p=p1, T=T1, v=vi, X0=X0, ceq=ceq)
    print("preshock:\n", preshock, '\n')

    postshock = normal_shock(preshock)
    print("postshock:\n", postshock, '\n')

    ionisation_fraction = postshock.X[spnames.index('e-')]*2
    print("Ionization fraction: {:g} %".format(ionisation_fraction*100))

    Nav = 6.02214076e+23
    Y = ceq.XtoY(postshock.X)
    eidx = spnames.index('e-')
    Me = ceq.M[eidx]
    rhoe = postshock.rho*Y[eidx]
    ne = rhoe*Nav/Me
    print("Electron Number Density: {:e}".format(ne))

    # Compare to CEA shock problem:
    # MOLE FRACTIONS
    output = """
         U2, M/SEC         474.76
         RHO, KG/CU M    8.9149-2
         P, BAR            2.3447
         T, K             6568.51
        *e-             3.5760-4
        *N              2.3939-1
        *N+             1.0167-5
        *NO             8.5898-3
        *NO+            3.2993-4
        *N2             4.2615-1
        *N2+            4.0200-6
        *O              3.2489-1
        *O+             1.4390-5
        *O-             1.0156-6
        *O2             2.6103-4
    """
    ceadata = {}
    for line in output.splitlines():
        line = line.strip()
        if len(line)==0: continue

        items = line.strip().split()
        k = items[0]
        v = items[-1]
        spname = k.strip().lstrip('*').rstrip(',')
        value = float(v.strip().replace('-','e-'))
        ceadata[spname] = value

    formatdiff = lambda a,b : "{:<9.1f}  {:<9.1f}  ({:<4.2f} %)".format(a, b, abs(a-b)/a*100)
    eformatdiff = lambda a,b : "{:<8.3e}  {:<8.3e}  ({:<4.2f} %)".format(a, b, abs(a-b)/a*100)
    print("Postshock Gas State:")
    print("           ceq        CEA        (% diff)")
    print("-------------------------------------------")
    print(" v (m/s)|  {}".format(formatdiff(postshock.v, ceadata['U2'])))
    print(" T (K)  |  {}".format(formatdiff(postshock.T, ceadata['T'])))
    print(" p (Pa) |  {}".format(formatdiff(postshock.p, ceadata['P']*1e5)))

    sp_ordered = [name for _,name in sorted(zip(postshock.X, spnames), reverse=True)]
    for spname in sp_ordered:
        idx = spnames.index(spname)
        print(" X {:3}  |  {}".format(spname, eformatdiff(postshock.X[idx], ceadata[spname])))

