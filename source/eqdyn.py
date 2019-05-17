"""
ESTCj-like functionality without cea

@author: Nick Gibbons
"""

from numpy import sum
from pyeq import *
from gasdynamics import IdealGas

Ru = 8.3144598

class EquilibriumGas(object):
    def __init__(self,rho=None,T=None,p=None,v=None,rhot=None,Tt=None,
                 pt=None,M=None,theta=0.0,spnames=None,massf0=None):
        
        self.rho =rho; self.T=T; self.p=p; self.v=v;
        self.rhot=rhot; self.Tt=Tt; self.pt=pt; self.M=M
        self.theta=theta; self.spnames = spnames; self.massf0 = massf0

        elements, nsp, nel, lewisdata, a, Ms = startup(species)
        self.elements = elements
        self.nsp      = nsp
        self.nel      = nel
        self.lewisdata= lewisdata
        self.a        = a
        self.Ms       = Ms
        self.lib      = load_ceq_library()

        # For starters assume that T, p are specified.
    
    @staticmethod 
    def Mmix_from_molef(X,MS):
        return sum(X*Ms)

    @staticmethod 
    def Mmix_from_massf(X,MS):
        return sum(Y/Ms)**-1

    @staticmethod
    def massf_from_molef(X,Ms):
        Mmix = sum(X*Ms)
        return X*Ms/Mmix

    @staticmethod
    def molef_from_massf(Y,Ms):
        Mmix = sum(Y/Ms)**-1
        return Y*Mmix/Ms

    # You could cache these maybe
    @property
    def R(self):
        return Ru/self.Mmix_from_massf(self.Y, self.Ms)

    @property
    def cp(self):
        return None

    @property
    def k(self):
        return None # This requires knowledge of T. What if we don't know it??? Iterate??? Shit. 

if __name__=='__main__':
    Air = EquilibriumGas(T=600.0, p=101.35, v=100.0, spnames=['N2','O2'], massf0=[0.75,0.25])
     
