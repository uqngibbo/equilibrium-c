/*
Library for equilibrium chemistry calculations, ported from ceq.c
    References:
        "Computer Program for Calculation of Complex Equilibrium Compositions and Applications"
        Nasa Reference Publication 1311, October 1995
        Sanford Gordon and Bonnie J. McBride

        "NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species"
        NASA/TP - 2002-211556, September 2002
        Bonnie J. McBride, Michael J. Zehe, and Sanford Gordon

@author: Nick Gibbons
*/

module ceq;
import std.stdio;
import std.math;

const double Ru=8.3144621;      // Universal Gas Constant (J/mol/K)
const double TRACELIMIT=1e-7;   // Trace species limiter (for ns/n)

@nogc double get_u(double T, double* X, int nsp, double* lewis, double* M){
    /*
    Compute thermal equilibrium u from known composition and primitives
    Inputs:
        T     : Temperature (K)
        X     : Composition [nsp]
        nsp   : number of species 
        lewis : Nasa Lewis Thermodynamic Database Data [nsp*3*9]
        M     : Molar Mass of each species (kg/mol) [nsp]
        verbose: print debugging information

    Output:
        u : internal energy per unit mass
    */
    int s;
    double Mmix, u, ns, U0_RTs;
    double* lp;

    Mmix = 0.0; for (s=0; s<nsp; s++) Mmix+=X[s]*M[s];
    
    u = 0.0;
    for (s=0; s<nsp; s++){
        ns = X[s]/Mmix;
        lp = lewis + 9*3*s;
        U0_RTs = compute_H0_RT(T, lp) - 1.0;
        u += ns*U0_RTs*Ru*T;
    }
    return u;
}

@nogc double get_cp(double T, double* X, int nsp, double* lewis, double* M){
    /*
    Compute thermal equilibrium cp from known composition and primitives
    Inputs:
        T     : Temperature (K)
        X     : Composition [nsp]
        nsp   : number of species 
        lewis : Nasa Lewis Thermodynamic Database Data [nsp*3*9]
        M     : Molar Mass of each species (kg/mol) [nsp]
        verbose: print debugging information

    Output:
        cp : specific heat at constant pressure per unit mass
    */
    int s;
    double Mmix, ns, cp,Cp0_Rs;
    double* lp;

    Mmix = 0.0; for (s=0; s<nsp; s++) Mmix+=X[s]*M[s];
    
    cp = 0.0;
    for (s=0; s<nsp; s++){
        ns = X[s]/Mmix;
        lp = lewis + 9*3*s;
        Cp0_Rs = compute_Cp0_R(T, lp);
        cp += ns*Cp0_Rs*Ru;
    }
    return cp;
}

// static confines this function to this module
//void pt_Assemble_Matrices(double* a,double* bi0,double* G0_RTs,double p,double* ns,double n,
//                          int nsp,int nel,double* A, double* B){
//    /*
//    Construct Iteration Matrix for reduced Newton Rhapson step, (eqn 2.24 and 2.26 from cea_I)
//    */
//    double lnns, lnn, lnp, akjaijnj, akjnjmuj, mus_RTj, bk;
//    double nss, nsmus;
//    int k,neq,i,j,s;
//    neq = nel+1;
//    lnn = log(n);
//    lnp = log(p/1e5); // Standard pressure for the tables is one BAR
//
//    // Equation 2.24: k-> equation index, i-> variable index
//    for (k=0; k<nel; k++){
//        bk = 0.0; for (s=0; s<nsp; s++) bk += a[k*nsp + s]*ns[s];
//        A[k*neq + 0] = bk;
//
//        for (i=0; i<nel; i++){
//            akjaijnj = 0.0;
//            for (j=0; j<nsp; j++){
//                akjaijnj += a[k*nsp+j]*a[i*nsp+j]*ns[j];
//            }
//            A[k*neq + i+1] = akjaijnj;
//        }
//        akjnjmuj = 0.0;
//        for (j=0; j<nsp; j++){
//            if (ns[j]==0.0) continue;
//            mus_RTj = G0_RTs[j] + log(ns[j]) - lnn + lnp;
//            akjnjmuj += a[k*nsp+j]*ns[j]*mus_RTj;
//
//        }
//        B[k] = bi0[k] - bk + akjnjmuj;
//
//        // Equation 2.26 - > (only the pii entries, we're highjacking k here to go across the last row)
//        A[nel*neq + k+1] = bk;
//    }
//    // Equation 2.26 - > (now the rest)
//    nss = 0.0;
//    nsmus = 0.0;
//    for (j=0; j<nsp; j++){
//        if (ns[j]==0.0) continue;
//        mus_RTj = G0_RTs[j] + log(ns[j]) - lnn + lnp; // I guess its okay to compute this again
//        nss += ns[j];
//        nsmus += ns[j]*mus_RTj;
//    }
//    A[nel*neq + 0]  = nss - n; // Yes, this is a problem. Shieeeeeet.
//    B[nel] = n - nss + nsmus;
//
//    //for (i=0; i<neq; i++){
//    //    for (j=0; j<neq; j++){
//    //        printf("%f ", A[i*neq+j]);
//    //    }
//    //    printf("| %f\n", B[i]);
//    //}
//    return;
//}
//
//void pt_species_corrections(double* S,double* a,double* G0_RTs,double p,double n,double* ns,
//                        int nsp, int nel, double* dlnns){
//    /*
//    Compute delta_log(ns) from the reduced iteration equations from 
//    equation 2.18m using the other deltas in S
//    Inputs:
//        S      : Corrections (pi1, pi2, pi3 ... dlog(n) [nel+1]
//        a      : elemental composition array [nel,nsp]
//        G0_RTs : Gibbs free energy of each species, divided by RT [nsp]
//        p      : pressure 
//        n      : total moles/mixture kg 
//        ns     : species moles/mixture kg [nsp]
//        nsp    : total number of species
//        nel    : total  number of elements 
//
//    Outputs:
//        dllns : change in log(ns) [nsp]
//    */
//    double dlnn,aispii,mu_RTs,lnn,lnp;
//    int s,i;
//    dlnn = S[0];
//    lnn = log(n);
//    lnp = log(p/1e5);
//
//    for (s=0; s<nsp; s++) {
//        if (ns[s]==0.0) { dlnns[s] = 0.0; continue;}
//        mu_RTs = G0_RTs[s] + log(ns[s]) - lnn + lnp;
//
//        aispii = 0.0;
//        for (i=0; i<nel; i++){
//            aispii += a[i*nsp+s]*S[i+1]; // S[i+1] = pi_i, the lagrange multiplier
//        }
//        dlnns[s] = -mu_RTs + dlnn + aispii;
//        //printf("dllns[%d] %f %f %f (%f) %f\n",s,-mu_RTs, dlnn, aispii, dlnns[s], lnn);
//    }
//    return; 
//}
//
//void pt_update_unknowns(double* S,double* dlnns,int nsp,double* ns,double* n){
//    /*
//    Add corrections to unknown values (ignoring lagrange multipliers)
//    Inputs:
//        S : vector of corrections from matrix step [nel+1]
//        dlnns : vector of species mole/mixture corrections [nsp]
//        nsp : number of species
//    Outputs:
//        ns : vector of species mole/mixtures [nsp]
//        n  : pointer to total moles/mixture (passed by reference!) [1]
//    */
//    int s;
//    double lnns,lnn,n_copy,lambda;
//
//    lnn = log(*n); // compute the log of the thing n is pointing to
//    lambda = fmin(1.0, 0.5*fabs(lnn)/fabs(S[0]));
//    n_copy = exp(lnn + lambda*S[0]); 
//    *n = n_copy;   // thing pointed to by n set to exp(lnn + S[0]);
//
//    for (s=0; s<nsp; s++){
//        if (ns[s]==0.0) continue;
//        lnns = log(ns[s]);
//
//        lambda = fmin(1.0, fabs(lnn)/fabs(dlnns[s]));
//        ns[s] = exp(lnns + lambda*dlnns[s]);
//
//        if (ns[s]/n_copy<TRACELIMIT) ns[s] = 0.0;
//    }
//
//    return;
//}
//
//// FIXME: This won't work because the calling function is @nogc and your working arrays are being
//// allocated in here. You'll need to put them in the gas model private attributes, possible also 
//// moving this function too.
//int solve_pt(double p,double T,double* X0,int nsp,int nel,double* lewis,double* M,double* a,double* X1,int verbose){
//    /*
//    Compute the equilibrium composition X1 at a fixed temperature and pressure
//    Inputs:
//        p     : Pressure (Pa)
//        T     : Temperature (K)
//        X0    : Intial Mole fractions [nsp]
//        nsp   : number of species 
//        nel   : number of elements 
//        lewis : Nasa Lewis Thermodynamic Database Data [nsp*3*9]
//        M     : Molar Mass of each species (kg/mol) [nsp]
//        a     : elemental composition array [nel,nsp]
//        verbose: print debugging information
//
//    Output:
//        X1 : Equilibrium Mole Fraction [nsp]  
//    */
//    double[] A_, B_, S_, G0_RTs_, ns_, bi0_, dlnns_; // Dynamic arrays
//    double* A, B, S, G0_RTs, ns, bi0, dlnns; // Pointers to preserve c mechanics
//    double* lp;
//    int neq,s,i,k,ntrace,errorcode;
//    double M0,n,M1,errorL2,thing;
//
//    const double tol=1e-6;
//    const int attempts=30;
//
//    errorcode=0;
//    neq= nel+1;
//    A_.length      = neq*neq; // Iteration Jacobian
//    B_.length      = neq;     // Iteration RHS
//    S_.length      = neq;     // Iteration unknown vector
//    G0_RTs_.length = nsp;     // Species Gibbs Free Energy
//    ns_.length     = nsp;     // Species moles/mixture mass
//    bi0_.length    = nel;     // starting composition coefficients
//    dlnns_.length  = nsp;     // starting composition coefficients
//
//    A      = A_.ptr;
//    B      = B_.ptr;
//    S      = S_.ptr;
//    G0_RTs = G0_RTs_.ptr;
//    ns     = ns_.ptr;
//    bi0    = bi0_.ptr;
//    dlnns  = dlnns_.ptr;
//
//    // Initialise Arrays and Iteration Guesses
//    M0 = 0.0;
//    for (s=0; s<nsp; s++) M0 += M[s]*X0[s];
//    for (i=0; i<nel; i++){
//        bi0[i] = 0.0;
//        for (s=0; s<nsp; s++){
//            bi0[i] += a[i*nsp + s]*X0[s]/M0;
//            }
//    }
//    n = 0.0;
//    for (s=0; s<nsp; s++) n += X0[s]/M0;
//    for (s=0; s<nsp; s++) ns[s] = n/nsp;
//    n*=1.1;
//
//    for (s=0; s<nsp; s++){
//        lp = lewis + 9*3*s;
//        G0_RTs[s] = compute_G0_RT(T, lp);
//    }
//
//    // Main Iteration Loop: 
//    for (k=0; k<=attempts; k++){
//        // 1: Perform an update of the equations
//        pt_Assemble_Matrices(a, bi0, G0_RTs, p, ns, n, nsp, nel, A, B);
//        solve_matrix(A, B, S, neq);
//        pt_species_corrections(S, a, G0_RTs, p, n, ns, nsp, nel, dlnns);
//        pt_update_unknowns(S, dlnns, nsp, ns, &n);
//
//        // Compute remaining error by checking species corrections
//        errorL2 = 0.0;
//        for (s=0; s<nsp; s++) errorL2 += dlnns[s]*dlnns[s];
//        errorL2 = pow(errorL2, 0.5);
//        if (verbose>0) printf("iter %d: %f %f %f %f  (%e) \n", k, n, ns[0], ns[1], ns[2], errorL2);
//
//        // Check for imminent singularity
//        ntrace=0;
//        for (s=0; s<nsp; s++) if (ns[s]==0.0) ntrace++;
//        if (ntrace==nsp-1) { // Pseudo convergence criteria, all the species but one are trace
//            for (s=0; s<nsp; s++) if (ns[s]!=0.0) i=s;
//            n = 1.0/M[i];
//            ns[i] = n;
//            errorL2 = 0.0;
//            if (verbose>0) printf("Pseudo convergence! Remaining species: (%d)\n", i);
//        }
//
//        // Exit loop of convergence is achieved
//        if (errorL2<tol) break;
//
//        if (isNaN(errorL2)) {
//            printf("Solver nan'd, exiting!\n");
//            errorcode=1;
//            break;
//        }
//
//        // Exit loop if too many attempts are undertaken
//        if (k==attempts) {
//            printf("Solver not converged, exiting!\n");
//            errorcode=1;
//            break;
//        }
//    }
//    
//    if ((verbose>0)&&(errorcode==0)) printf("Converged in %d iter, error: %e\n", k, errorL2);
//    // Compute output composition
//    M1 = 1.0/n;
//    for (s=0; s<nsp; s++) X1[s] = M1*ns[s];
//
//    // Don't need to clear memory manually
//    return errorcode;
//}

//void rhou_Assemble_Matrices(double* a,double* bi0, double rho0,double u0,double T,double* ns,int nsp,
//                              int nel,double* A, double* B, double* G0_RTs, double* U0_RTs,
//                              double* Cv0_Rs, double* lewis){
//    /*
//    Construct Iteration Matrix for reduced Newton Rhapson rhou step, (eqn 2.45 and 2.47 from cea_I)
//    */
//    double mus_RTj, bk, u, akjaijnj, akjnjmuj, akjnjUj, njCvj, njUj2, njUjmuj, aijnjUj;
//    int k,neq,i,j,s;
//    double *lp;
//    neq = nel+1;
//
//    u = 0.0;
//    for (s=0; s<nsp; s++){
//        lp = lewis + 9*3*s;
//        G0_RTs[s] = compute_G0_RT(T, lp);
//        U0_RTs[s] = compute_H0_RT(T, lp) - 1.0;
//        Cv0_Rs[s] = compute_Cp0_R(T, lp) - 1.0;
//        u += ns[s]*U0_RTs[s]*Ru*T;
//    }
//
//    // Equation 2.45: k-> equation index, i-> variable index
//    for (k=0; k<nel; k++){
//        bk = 0.0; for (s=0; s<nsp; s++) bk += a[k*nsp + s]*ns[s];
//
//        for (i=0; i<nel; i++){
//            akjaijnj = 0.0;
//            for (j=0; j<nsp; j++){
//                akjaijnj += a[k*nsp+j]*a[i*nsp+j]*ns[j];
//            }
//            A[k*neq + i+1] = akjaijnj; // Lagrange multiplier matrix entry
//        }
//        akjnjmuj = 0.0;
//        akjnjUj = 0.0;
//        for (j=0; j<nsp; j++){
//            mus_RTj = G0_RTs[j] + log(rho0*ns[j]*Ru*T/1e5);
//            akjnjmuj += a[k*nsp+j]*ns[j]*mus_RTj;
//            akjnjUj  += a[k*nsp+j]*ns[j]*U0_RTs[j];
//        }
//        A[k*neq + 0] = akjnjUj; // Temperature matrix entry
//        B[k] = bi0[k] - bk + akjnjmuj; // RHS of kth equation 2.45
//
//    }
//    // Equation 2.47
//    njCvj  = 0.0;
//    njUj2  = 0.0;
//    njUjmuj= 0.0;
//
//    for (j=0; j<nsp; j++){
//        mus_RTj = G0_RTs[j] + log(rho0*ns[j]*Ru*T/1e5);
//        njCvj += ns[j]*Cv0_Rs[j];
//        njUj2 += ns[j]*U0_RTs[j]*U0_RTs[j];
//        njUjmuj += ns[j]*U0_RTs[j]*mus_RTj;
//    }
//    A[nel*neq + 0]  = njCvj + njUj2;  // Temperature matrix entry
//    B[nel] = (u0 - u)/Ru/T + njUjmuj; // RHS 
//
//    for (i=0; i<nel; i++){
//        aijnjUj = 0.0;
//        for (j=0; j<nsp; j++){
//            aijnjUj += a[i*nsp + j]*ns[j]*U0_RTs[j];
//        }
//        A[nel*neq + 1+i] = aijnjUj; // Lagrange Multipliers       
//    }
//
//    //for (i=0; i<neq; i++){
//    //    for (j=0; j<neq; j++){
//    //        printf("%f ", A[i*neq+j]);
//    //    }
//    //    printf("| %f\n", B[i]);
//    //}
//    return;
//}
//
//void rhou_species_corrections(double* S,double* a,double* G0_RTs,double* U0_RTs,double rho0,double T,
//                         double* ns, int nsp, int nel, double* dlnns){
//    /*
//    Compute delta_log(ns) from the reduced iteration equations from 
//    equation 2.18m using the other deltas in S
//    Inputs:
//        S      : Corrections (pi1, pi2, pi3 ... dlog(n) [nel+1]
//        a      : elemental composition array [nel,nsp]
//        G0_RTs : Gibbs free energy of each species, divided by RT [nsp]
//        U0_RTs : Internal energy of each species, divided by RT [nsp]
//        rho0   : goal density (kg/m3)
//        T      : current temperature guess (K)
//        ns     : species moles/mixture kg [nsp]
//        nsp    : total number of species
//        nel    : total  number of elements 
//
//    Outputs:
//        dllns : change in log(ns) [nsp]
//    */
//    double dlnT,mus_RTs,aispii;
//    int s,i;
//    dlnT = S[0];
//
//    // These should (should?) have been set during the assemble matrix step
//    //for (s=0; s<nsp; s++){
//    //    lp = lewis + 9*3*s;
//    //    G0_RTs[s] = compute_G0_RT(T, lp);
//    //    U0_RTs[s] = compute_H0_RT(T, lp) - 1.0;
//    //}
//    
//    for (s=0; s<nsp; s++) {
//        mus_RTs = G0_RTs[s] + log(rho0*ns[s]*Ru*T/1e5);
//
//        aispii = 0.0;
//        for (i=0; i<nel; i++){
//            aispii += a[i*nsp+s]*S[i+1]; // S[i+1] = pi_i, the lagrange multiplier
//        }
//        dlnns[s] = -mus_RTs + U0_RTs[s]*dlnT + aispii;
//    }
//    return; 
//}
//
//void rhou_update_unknowns(double* S,double* dlnns,int nsp,double* ns,double* T){
//    /*
//    Add corrections to unknown values (ignoring lagrange multipliers)
//    Inputs:
//        S : vector of corrections from matrix step [nel+1]
//        dlnns : vector of species mole/mixture corrections [nsp]
//        nsp : number of species
//    Outputs:
//        ns : vector of species mole/mixtures [nsp]
//        T  : pointer to temperature guess (passed by reference!) [1]
//    */
//    int s;
//    double lnns,lnT,n,lnn,lambda;
//    lnT = log(*T); // compute the log of the thing T is pointing to
//    lambda = fmin(1.0, 0.5*fabs(lnT)/fabs(S[0]));
//    *T = exp(lnT + lambda*S[0]); // thing pointed to by T set to exp(lnT + S[0]);
//
//    n = 0.0; for (s=0; s<nsp; s++) n+=ns[s]; lnn=log(n);
//
//    for (s=0; s<nsp; s++){
//        lnns = log(ns[s]);
//        lambda = fmin(1.0, fabs(lnn)/fabs(dlnns[s]));
//        ns[s] = exp(lnns + lambda*dlnns[s]);
//
//        if (ns[s]/n<TRACELIMIT) ns[s] = 0.0;
//    }
//    return;
//}
//
//double temperature_guess(int nsp, double u, double M0, double* X0, double* lewis){
//    /*
//    Guess a first iteration temperature assuming constant Cv from 298 K
//    Inputs:
//        nsp   : Number of species
//        u     : Target mixture internal energy (J/kg)
//        M0    : Initial composition molecular weight (kg/mol)
//        X0    : Intiial composition mole fractions [nsp]
//        lewis : Nasa Lewis Thermodynamic Database Data [nsp*3*9]
//
//    Output:
//        T : Temperature Guess (K)
//    */
//    int s;
//    double* lp;
//    double uf,cv,T,ufs,cvs,Cps298,Hfs298,ns0;
//
//    uf = 0.0;
//    cv = 0.0;
//    for (s=0; s<nsp; s++){
//        lp = lewis + 9*3*s;
//        Cps298 = compute_Cp0_R(298.15, lp)*Ru;
//        Hfs298 = compute_H0_RT(298.15, lp)*Ru*298.15;
//
//        ns0 = X0[s]/M0;
//        ufs= ns0*(Cps298*298.15 - Hfs298);
//        cvs= ns0*(Cps298 - Ru);
//
//        uf += ufs;
//        cv += cvs;
//    }
//    T = (u + uf)/cv;
//    return T;
//}
//
//int solve_rhou(double rho,double u,double* X0,int nsp,int nel,double* lewis,double* M,double* a,
//               double* X1, double* Teq, int verbose){
//    /*
//    Compute the equilibrium composition X1 at a fixed volume and internal energy 
//    Inputs:
//        rho   : target Density (kg/m3)
//        u     : target internal energy (J/kg)
//        X0    : Intial Mole fractions [nsp]
//        nsp   : number of species 
//        nel   : number of elements 
//        lewis : Nasa Lewis Thermodynamic Database Data [nsp*3*9]
//        M     : Molar Mass of each species (kg/mol) [nsp]
//        a     : elemental composition array [nel,nsp]
//        verbose: print debugging information
//
//    Output:
//        X1  : Equilibrium Mole Fraction [nsp]  
//        Teq: Equilibrium Temperature 
//    */
//    double[] A_, B_, S_, G0_RTs_, U0_RTs_, Cv0_Rs_, ns_, bi0_, dlnns_; // Dynamic arrays
//    double* A, B, S, G0_RTs, U0_RTs, Cv0_Rs, ns, bi0, dlnns; // Pointers to preserve c mechanics
//    int neq,s,i,k,errorcode;
//    double M0,n,M1,errorL2,thing,T;
//
//    const double tol=1e-6;
//    const int attempts=10;
//
//    errorcode=0;
//    neq= nel+1;
//    A_.length      = neq*neq; // Iteration Jacobian
//    B_.length      = neq;     // Iteration RHS
//    S_.length      = neq;     // Iteration unknown vector
//    G0_RTs_.length = nsp;     // Species Gibbs Free Energy
//    U0_RTs_.length = nsp;     // Species Internal Energy
//    Cv0_Rs_.length = nsp;     // Species Specific Heat @ Const_. volume 
//    ns_.length     = nsp;     // Species moles/mixture mass
//    bi0_.length    = nel;     // starting composition coefficients
//    dlnns_.length  = nsp;     // starting composition coefficients
//
//    A      = A_.ptr;      
//    B      = B_.ptr;     
//    S      = S_.ptr;      
//    G0_RTs = G0_RTs_.ptr; 
//    U0_RTs = U0_RTs_.ptr; 
//    Cv0_Rs = Cv0_Rs_.ptr; 
//    ns     = ns_.ptr;     
//    bi0    = bi0_.ptr;    
//    dlnns  = dlnns_.ptr;
//
//    // Initialise Arrays and Iteration Guesses
//    M0 = 0.0;
//    for (s=0; s<nsp; s++) M0 += M[s]*X0[s];
//    for (i=0; i<nel; i++){
//        bi0[i] = 0.0;
//        for (s=0; s<nsp; s++){
//            bi0[i] += a[i*nsp + s]*X0[s]/M0;
//            }
//    }
//    n = 0.0;
//    for (s=0; s<nsp; s++) n += X0[s]/M0;
//    for (s=0; s<nsp; s++) ns[s] = n/nsp;
//
//    T = temperature_guess(nsp, u, M0, X0, lewis);
//    if (verbose>0) printf("Guess T: %f\n", T);
//
//    // Begin Iterations
//    for (k=0; k<attempts; k++){
//        rhou_Assemble_Matrices(a,bi0,rho,u,T,ns,nsp,nel,A,B,G0_RTs,U0_RTs,Cv0_Rs,lewis);
//        solve_matrix(A, B, S, neq);
//        rhou_species_corrections(S,a,G0_RTs,U0_RTs,rho,T,ns,nsp,nel,dlnns);
//        rhou_update_unknowns(S, dlnns, nsp, ns, &T);
//
//        errorL2 = 0.0;
//        for (s=0; s<nsp; s++) errorL2 += dlnns[s]*dlnns[s];
//        errorL2 = pow(errorL2, 0.5);
//        if (verbose>0) printf("iter %d: %f %f %f %f  (%e) \n", k, T, ns[0], ns[1], ns[2], errorL2);
//        if (errorL2<tol) break;
//
//        if (k>=attempts) {
//            printf("Solver not converged, exiting!\n");
//            errorcode=1;
//            break;
//        }
//    }
//    
//    if ((verbose>0)&&(errorcode==0)) printf("Converged in %d iter, error: %e\n", k, errorL2);
//    // Compute output composition
//    n = 0.0;
//    for (s=0; s<nsp; s++) n += ns[s];
//    M1 = 1.0/n;
//    for (s=0; s<nsp; s++) X1[s] = M1*ns[s];
//    *Teq = T;
//
//    // Don't need to clear memory manually
//    return errorcode;
//}


@nogc double compute_Cp0_R(double Tin, double* lewis){
    /*
    Molar Specific Heat at constant Pressure at 1 BAR divided by Ru*T for a species using
    the nasa lewis therodynamic data
     Tin   : Temperature (double)
     lewis : thermo data for current species (double)[3*9] 
    */
    double* lp;
    double Cp0_R,T;
    int iT;

    T = fmax(fmin(Tin,20000.0),200.0);
    iT = T <= 1000.0 ? 0 : (T<=6000 ? 1 : 2);
    lp = lewis + iT*9;

    Cp0_R =  lp[0]/T/T
           + lp[1]/T
           + lp[2]
           + lp[3]*T
           + lp[4]*T*T
           + lp[5]*T*T*T
           + lp[6]*T*T*T*T;
      return Cp0_R;
}

@nogc double compute_H0_RT(double Tin, double* lewis){
    /*
    Molar enthalpy at 1 BAR divided by Ru*T for a species using the nasa lewis therodynamic data
     Tin   : Temperature (double)
     lewis : thermo data for current species (double)[3*9] 
    */
    double* lp;
    double H0_RT,T;
    int iT;

    T = fmax(fmin(Tin,20000.0),200.0);
    iT = T <= 1000.0 ? 0 : (T<=6000 ? 1 : 2);
    lp = lewis + iT*9;

    H0_RT =    -1.0*lp[0]/T/T
           + log(T)*lp[1]/T
           +    1.0*lp[2]
           +    0.5*lp[3]*T
           +  1/3.0*lp[4]*T*T
           +   0.25*lp[5]*T*T*T
           +    0.2*lp[6]*T*T*T*T + lp[7]/T;
      return H0_RT;
}

@nogc double compute_S0_R(double Tin, double* lewis){
    /*
    Molar entropy at 1 BAR divided by Ru for a species using the nasa lewis therodynamic data
     Tin   : Temperature (double)
     lewis : thermo data for current species (double)[3*9] 
    */
    double* lp;
    double S0_R,T;
    int iT;

    T = fmax(fmin(Tin,20000.0),200.0);
    iT = T <= 1000.0 ? 0 : (T<=6000 ? 1 : 2);
    lp = lewis + iT*9;

    S0_R =     -0.5*lp[0]/T/T
           +   -1.0*lp[1]/T
           + log(T)*lp[2]
           +    1.0*lp[3]*T
           +    0.5*lp[4]*T*T
           +  1/3.0*lp[5]*T*T*T
           +   0.25*lp[6]*T*T*T*T + lp[8];
    return S0_R;
}

@nogc double compute_G0_RT(double T, double* lewis){
    /*
    Compute the Molar entropy divided by Ru*T for a species using the nasa lewis therodynamic data
     T     : Temperature (double)
     lewis : thermo data for current species (double)[3*9] 
    */
    return compute_H0_RT(T, lewis) - compute_S0_R(T, lewis);
}

@nogc void solve_matrix(double* A, double* B, double *X, int N){
    /*
    Bare bones Gaussian matrix solver for A*X = B 
       A : pointer to NxN matrix of doubles
       B : pointer to N vector of doubles
       X : pointer to N solution vector of doubles
       N : problem size (int)
    */
    int i,j,k,p;
    double pval,a,b;
    k=0; // working diagonal for row reduction

    while (k<N) {

    // Find the partial-pivot index
    pval = -1.0;
    p = k;
    for (i=k; i<N; i++) {
        a = fabs(A[i*N+k]);
        if (a>pval){
           pval=a;
           p=i;
           }
    }

    // Row swap
    for (j=0; j<N; j++) {
        a = A[k*N+j];
        A[k*N+j] = A[p*N+j];
        A[p*N+j] = a;
    }
    a = B[k];
    B[k] = B[p];
    B[p] = a;

    // Multiply Pivot Row
    a = 1.0/A[k*N+k];
    for (j=0; j<N; j++) {
        A[k*N+j]*=a;
    }
    B[k]*=a;

    // Add multiples of pivot row to descending rows
    for (i=k+1; i<N; i++) { // should skip if k==N
        a = A[i*N+k];
        for (j=0; j<N; j++){
            A[i*N+j] -= a*A[k*N+j];
        }
        B[i] -= a*B[k];
    }
    k++;

    } // end while loop for row reduction

    //for (i=0; i<N; i++) {
    //    printf("[");
    //    for (j=0; j<N; j++){
    //        printf("%f ", A[i*N+j]);
    //    }
    //    printf("]\n");
    //}
    //printf("\n");
    // Now compute X using back substitution
    for (k=N-1; k>=0; k--){
        a = 0.0;
        for (j=k+1; j<N; j++) a+= A[k*N+j]*X[j]; // Should skip when k=N
        X[k] = B[k] - a;
    }
    return;
}


/*
int main(){
    printf("Called main in ceq.c!\n");
    return 0;
}
*/
