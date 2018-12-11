"""
Optimiser based approach to chemical equilibrium
 - Solve pt problem for a simple set of species

@author: Nick Gibbons
"""
from lewis_thermo import get_species
from numpy import array, zeros, exp, log
from numpy.linalg import solve

def pteq(spnames, p, T, Xs0,tol=1e-6):
    species = [get_species(sp) for sp in spnames]
    elements = set()
    for s in species:
        for k in s.atoms.keys(): elements.add(k)
    elements = list(elements)

    nsp = len(species)
    nel = len(elements)

    a = zeros((len(elements),len(species)))
    for i,s in enumerate(species):
        for k,v in s.atoms.items():
            j = elements.index(k)
            a[j,i] = v

    # Initial Guesses
    M0 = sum(X*sp.M for X,sp in zip(Xs0,species))
    ns0 = Xs0/M0
    bi0 = (a*ns0).sum(axis=1)

    ns= ns0*0.0 + ns0.sum()/ns0.size
    n  = 1.1*ns.sum() 
    pii = zeros(bi0.shape)

    Ru = 8.3144598 # uses gram-moles, or regular moles... to match M in lewis_thermo 
    G0_RTs = [s.H0onRT(T) - s.S0onR(T) for s in species]

    for k in range(20):
        A,B = Assemble_Matrices(a, bi0, G0_RTs, p, ns, n)
        corrections = solve(A,B)
        error = corrections[:nsp+1]
        errorL2 = ((error.dot(error)).sum())**0.5
        ns, n, pii = update_unknowns(corrections, ns, n)

        if errorL2<tol: break
    else:
        print("Solver Not Converging!")
        print(ns, n, pii,corrections,errorL2)
        raise Exception

    return ns, n

if __name__=='__main__':
    spnames = ['CO2', 'CO', 'O2']
    T = 2500.0
    p = 0.1*101.35e3
    Xs0 = array([1.0, 0.0, 0.0])

    species = [get_species(sp) for sp in spnames]
    Ru = 8.3144598 # uses gram-moles, or regular moles... to match M in lewis_thermo 
    Xst = array([0.66108962603325838,0.22594024931116111,0.11297012465558055]) # target composition from co2_new.py 
    Mt = sum(X*sp.M for X,sp in zip(Xst,species))
    Rt = Ru/Mt
    rhot = p/Rt/T
    nt = 1/Mt
    cst = Xst*rhot/Mt
    nst = Xst/Mt # also cs/rhot 
    nt2 = nst.sum()
    print("Mt", Mt)
    print("rhot", rhot)
    print("cst"  , cst)
    print("nst"  , nst)
    print("nt"  , nt)
    print("nt2"  , nt2)

    ns, n = pteq(spnames, p, T, Xs0)
    M1 = 1.0/n
    Xs1 = M1*ns
    print("result", Xs1)
    print("target", Xst)
