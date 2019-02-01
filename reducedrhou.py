"""
Optimiser based approach to chemical equilibrium
 - Solve rho u problem for a simple set of species

 References:
  - nasacea_I

@author: Nick Gibbons
"""
from lewis_thermo import get_species
from numpy import array, zeros, exp, log
from numpy.linalg import solve

def Assemble_Matrices(a, bi0, rho0, u0, T, ns):
    """
    Assemble Jacobian Matrix and RHS for Newton-Rhapson iteration
    Inputs:
     - a : atomic element count (nel,nsp)
     - bi0 : starting atomic element count sum (nel)
     - rho0: target density 
     - u0  : target internal energy per unit mass
     - T   : Current Guess for the temperature
     - ns  : Array of current guesses for species n's (nsp)

    Output: 
     - A : Jacobian Matrix (nel+1, nel+1)
     - B : RHS Matrix (nel+1)
    """
    nsp= len(ns)
    nel= len(bi0)
    neq = nel + 1
    A = zeros((neq,neq))
    B = zeros((neq))

    G0_RTs = array([s.H0onRT(T) - s.S0onR(T) for s in species])
    U0_RTs = array([s.H0onRT(T) - 1.0 for s in species])
    Cv_Rs = array([s.Cp0onR(T) - 1.0 for s in species])

    Ru = 8.3144598 # uses gram-moles, or regular moles... to match M in lewis_library
    #n = ns.sum()
    #p = rho0*n*Ru*T
    #lnn = log(n)
    #lnns = log(ns)
    #lnp = log(p/1e5) # Because the standard pressure was one bar
    mu0s_RT = G0_RTs# - 1.0/rho0/Ru/T # ???? FIXME: What is this?
    mus_RT = mu0s_RT + log(rho0*ns*Ru*T/1e5) # Helmholtz chemical potential
    b = (a*ns).sum(axis=1)
    u = (ns*U0_RTs*Ru*T).sum()

    # Equation 2.45 
    for k in range(nel):
        for i in range(nel):
            akjaijnj = 0.0
            for j in range(nsp):
                akjaijnj+= a[k,j]*a[i,j]*ns[j]
            A[k,1+i] = akjaijnj # Lagrange multipliers

        akjnjmuj = 0.0
        akjnjUj = 0.0
        for j in range(nsp):
            akjnjmuj += a[k,j]*ns[j]*mus_RT[j]
            akjnjUj += a[k,j]*ns[j]*U0_RTs[j]

        A[k,0] = akjnjUj # Temperature correction
        B[k] = bi0[k] - b[k] + akjnjmuj 

    # Equation (2.47)
    njCvj = 0.0
    njUj2 = 0.0
    njUjmuj = 0.0
    for j in range(nsp):
        njCvj += ns[j]*Cv_Rs[j]
        njUj2 += ns[j]*U0_RTs[j]*U0_RTs[j]
        njUjmuj += ns[j]*U0_RTs[j]*mus_RT[j] 
    A[nel,0] = njCvj + njUj2 # Temperature correction

    for i in range(nel):
        aijnjUj = 0.0
        for j in range(nsp):
            aijnjUj += a[i,j]*ns[j]*U0_RTs[j]
        A[nel,1+i] = aijnjUj  # Lagrange multipliers

    B[nel] = (u0 - u)/Ru/T + njUjmuj 
    #print(A,B)

    return A,B

def get_dlnns(corrections, a, rho0, ns, T):
    """ Compute dln(ns) from equation 2.40 and the reduced system corrections """
    nel,nsp = a.shape 
    dlnT = corrections[0]
    pii = corrections[1:]

    G0_RTs = array([s.H0onRT(T) - s.S0onR(T) for s in species])
    U0_RTs = array([s.H0onRT(T) - 1.0 for s in species])
    #n = ns.sum()
    #p = rho0*n*Ru*T
    #lnn = log(n)
    #lnns = log(ns)
    #lnp = log(p/1e5) # Because the standard pressure was one bar
    mu0s_RT = G0_RTs# - 1.0/rho0/Ru/T # ???? which density do I use here?
    mus_RT = mu0s_RT + log(rho0*ns*Ru*T/1e5) # Helmholtz chemical potential

    aijpii = zeros(nsp) 
    for i in range(nel): aijpii += a[i,:]*pii[i]

    dlnns = -mus_RT + U0_RTs*dlnT + aijpii
    return dlnns

def update_unknowns(corrections, dlnns, ns, T):
    # unpack corrections and rescale
    dlnT = corrections[0]
    pii = corrections[1:]

    ns_new = exp(log(ns) + dlnns)
    T_new = exp(log(T) + dlnT)
    return T_new, ns_new, pii

def rhoueq_reduced(spnames, rho, u, Xs0,tol=1e-6):
    species = [get_species(sp) for sp in spnames]
    elements = set()
    for s in species:
        for k in s.atoms.keys(): elements.add(k)
    elements = list(elements)
    elements.sort()
    print(elements)

    nsp = len(species)
    nel = len(elements)

    a = zeros((len(elements),len(species)))
    for i,s in enumerate(species):
        for k,v in s.atoms.items():
            j = elements.index(k)
            a[j,i] = v

    M0 = sum(X*sp.M for X,sp in zip(Xs0,species))
    ns0 = Xs0/M0
    bi0 = (a*ns0).sum(axis=1)

    # Initial Guesses (change this to have the option to reuse Xs0)
    Ru = 8.3144598 # uses gram-moles, or regular moles... to match M in lewis_library
    pii = zeros(bi0.shape)
    ns= ns0*0.0 + ns0.sum()/ns0.size

    Cps298 = array([sp.Cp0onR(298.15)*Ru for sp in species])
    Hf298 = array([sp.H0onRT(298.15)*Ru*298.15 for sp in species])
    uf = ns0*(Cps298*298.15 - Hf298) # Use ns or ns0???
    cv = ns0*(Cps298 - Ru)
    T = (u + uf.sum())/(cv.sum()) # Guess using constant Cp at 298.15
    print("Guess T:", T)

    for k in range(20):
        A,B = Assemble_Matrices(a, bi0, rho, u, T, ns)
        corrections = solve(A,B)
        print("corrections:", corrections)
        dlnns = get_dlnns(corrections, a, rho, ns, T)
        error = dlnns
        errorL2 = ((error.dot(error)).sum())**0.5
        print("dlnns:", dlnns, errorL2)
        T, ns, pii = update_unknowns(corrections, dlnns, ns, T)
        print("Iter:", T, ns)
        print(' ')

        if errorL2<tol: break
    else:
        print("Solver Not Converging!")
        print(ns, n, pii,corrections,errorL2)
        raise Exception

    n = ns.sum()
    return ns, n


if __name__=='__main__':
    spnames = ['CO2', 'CO', 'O2']
    T = 2500.0
    p = 0.1*101.35e3
    Xs0 = array([1.0, 0.0, 0.0])
    Xst = array([0.66108962603325838,0.22594024931116111,0.11297012465558055]) # target composition from co2_new.py 

    species = [get_species(sp) for sp in spnames]
    Ru = 8.3144598 # uses gram-moles, or regular moles... to match M in lewis_thermo 
    Mt = sum(X*sp.M for X,sp in zip(Xst,species))
    Rt = Ru/Mt
    rhot = p/Rt/T
    nt = 1/Mt
    cst = Xst*rhot/Mt
    nst = Xst/Mt # also cs/rhot 
    nt2 = nst.sum()

    U0st= [(sp.H0onRT(T) - 1.0)*Ru*T for sp in species]
    ut = sum(njt*U0jt for njt,U0jt in zip(nst,U0st))

    print("Mt", Mt)
    print("rhot", rhot)
    print("cst" , cst)
    print("nst" , nst)
    print("nt"  , nt)
    print("nt2" , nt2)
    print("ut"  , ut)
    print("Ut"  , U0st)
    print("Tt"  , T)

    ns, n = rhoueq_reduced(spnames, rhot, ut, Xs0)
    M1 = 1.0/n
    Xs1 = M1*ns
    print("result", Xs1)
    print("target", Xst)
