/*
C library for equilibrium chemistry calculations
    References:
        "Computer Program for Calculation of Complex Equilibrium Compositions and Applications"
        Nasa Reference Publication 1311, October 1995
        Sanford Gordon and Bonnie J. McBride

        "NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species"
        NASA/TP - 2002-211556, September 2002
        Bonnie J. McBride, Michael J. Zehe, and Sanford Gordon

@author: Nick Gibbons
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "thermo.h"
#include "linalg.h"
#include "rhou.h"

static double constraint_errors(double* a,double* bi0,double* ns,int nsp,int nel,double* dlnns,int verbose){
    /*
    Unified computation of current error to determine whether to break loop
    */ 
    int s,i;
    double bi,errorrms,error;

    // Compute the change in current variables (note this is the unlimited dlnns! not the real change)
    errorrms = 0.0;

    for (s=0; s<nsp; s++) errorrms += dlnns[s]*dlnns[s];

    for (i=0; i<nel; i++) {
        bi = 0.0;
        for (s=0; s<nsp; s++) bi += a[i*nsp + s]*ns[s];
        error = bi - bi0[i];
        errorrms += error*error;
    }
    errorrms /= (nsp+nel);
    errorrms = sqrt(errorrms);
    return errorrms;
}

static void Assemble_Matrices(double* a,double* bi0, double rho0,double u0,double T,double* ns,int nsp,
                              int nel,double* A, double* B, double* G0_RTs, double* U0_RTs,
                              double* Cv0_Rs, double* lewis){
    /*
    Construct Iteration Matrix for reduced Newton Rhapson rhou step, (eqn 2.45 and 2.47 from cea_I)
    */
    double mus_RTj, bk, u, akjaijnj, akjnjmuj, akjnjUj, njCvj, njUj2, njUjmuj, aijnjUj;
    int k,neq,i,j,s;
    double *lp;
    neq = nel+1;

    u = 0.0;
    for (s=0; s<nsp; s++){
        lp = lewis + 9*3*s;
        G0_RTs[s] = compute_G0_RT(T, lp);
        U0_RTs[s] = compute_H0_RT(T, lp) - 1.0;
        Cv0_Rs[s] = compute_Cp0_R(T, lp) - 1.0;
        u += ns[s]*U0_RTs[s]*Ru*T;
    }

    // Equation 2.45: k-> equation index, i-> variable index
    for (k=0; k<nel; k++){
        bk = 0.0; for (s=0; s<nsp; s++) bk += a[k*nsp + s]*ns[s];

        for (i=0; i<nel; i++){
            akjaijnj = 0.0;
            for (j=0; j<nsp; j++){
                akjaijnj += a[k*nsp+j]*a[i*nsp+j]*ns[j];
            }
            A[k*neq + i+1] = akjaijnj; // Lagrange multiplier matrix entry
        }
        akjnjmuj = 0.0;
        akjnjUj = 0.0;
        for (j=0; j<nsp; j++){
            mus_RTj = G0_RTs[j] + log(rho0*ns[j]*Ru*T/1e5);
            akjnjmuj += a[k*nsp+j]*ns[j]*mus_RTj;
            akjnjUj  += a[k*nsp+j]*ns[j]*U0_RTs[j];
        }
        A[k*neq + 0] = akjnjUj; // Temperature matrix entry
        B[k] = bi0[k] - bk + akjnjmuj; // RHS of kth equation 2.45

    }
    // Equation 2.47
    njCvj  = 0.0;
    njUj2  = 0.0;
    njUjmuj= 0.0;

    for (j=0; j<nsp; j++){
        mus_RTj = G0_RTs[j] + log(rho0*ns[j]*Ru*T/1e5);
        njCvj += ns[j]*Cv0_Rs[j];
        njUj2 += ns[j]*U0_RTs[j]*U0_RTs[j];
        njUjmuj += ns[j]*U0_RTs[j]*mus_RTj;
    }
    A[nel*neq + 0]  = njCvj + njUj2;  // Temperature matrix entry
    B[nel] = (u0 - u)/Ru/T + njUjmuj; // RHS 

    for (i=0; i<nel; i++){
        aijnjUj = 0.0;
        for (j=0; j<nsp; j++){
            aijnjUj += a[i*nsp + j]*ns[j]*U0_RTs[j];
        }
        A[nel*neq + 1+i] = aijnjUj; // Lagrange Multipliers       
    }

    //for (i=0; i<neq; i++){
    //    for (j=0; j<neq; j++){
    //        printf("%f ", A[i*neq+j]);
    //    }
    //    printf("| %f\n", B[i]);
    //}
    return;
}

static void species_corrections(double* S,double* a,double* G0_RTs,double* U0_RTs,double rho0,double T,
                         double* ns, int nsp, int nel, double* dlnns, int verbose){
    /*
    Compute delta_log(ns) from the reduced iteration equations from 
    equation 2.18m using the other deltas in S
    Inputs:
        S      : Corrections (pi1, pi2, pi3 ... dlog(n) [nel+1]
        a      : elemental composition array [nel,nsp]
        G0_RTs : Gibbs free energy of each species, divided by RT [nsp]
        U0_RTs : Internal energy of each species, divided by RT [nsp]
        rho0   : goal density (kg/m3)
        T      : current temperature guess (K)
        ns     : species moles/mixture kg [nsp]
        nsp    : total number of species
        nel    : total  number of elements 

    Outputs:
        dllns : change in log(ns) [nsp]
    */
    double dlnT,mus_RTs,aispii;
    int s,i;
    dlnT = S[0];

    // These should (should?) have been set during the assemble matrix step
    //for (s=0; s<nsp; s++){
    //    lp = lewis + 9*3*s;
    //    G0_RTs[s] = compute_G0_RT(T, lp);
    //    U0_RTs[s] = compute_H0_RT(T, lp) - 1.0;
    //}
    
    for (s=0; s<nsp; s++) {
        mus_RTs = G0_RTs[s] + log(rho0*ns[s]*Ru*T/1e5);

        aispii = 0.0;
        for (i=0; i<nel; i++){
            aispii += a[i*nsp+s]*S[i+1]; // S[i+1] = pi_i, the lagrange multiplier
        }
        dlnns[s] = -mus_RTs + U0_RTs[s]*dlnT + aispii;
    }
    return; 
}

static void handle_singularity(double* S,double* a,double* G0_RTs,double* U0_RTs,double rho, double T,
                         double* ns, double n, int nsp, int nel, double* dlnns, int verbose){
    /*
    Compute delta_log(ns) from the reduced iteration equations from 
    equation 2.18m using the other deltas in S
    Inputs:
        S      : Corrections (pi1, pi2, pi3 ... dlog(n) [nel+1]
        a      : elemental composition array [nel,nsp]
        G0_RTs : Gibbs free energy of each species, divided by RT [nsp]
        U0_RTs : Internal energy of each species, divided by RT [nsp]
        rho0   : goal density (kg/m3)
        T      : current temperature guess (K)
        ns     : species moles/mixture kg [nsp]
        n      : mixture moles/mixture kg
        nsp    : total number of species
        nel    : total  number of elements 
        dllns : change in log(ns) [nsp]
    */
    const double RESET=4.0;
    int s;

    if (verbose>0) printf("    ### Singular Matrix!: Unlocking to %f\n",log(RESET*n*TRACELIMIT));

    for (s=0; s<nsp; s++){
        if (ns[s]!=0.0) continue;  // Ignore non trace species

        ns[s] = fmax(RESET*n*TRACELIMIT, ns[s]); // Reset trace trace species to a small but finite number
        species_corrections(S,a,G0_RTs,U0_RTs,rho,T,ns,nsp,nel,dlnns,0);
        if (dlnns[s]<0.0) ns[s] = 0.0; // Re-zero any species with negative predicted dlnns
        if (verbose>0) printf("   faux dlnns: %f changed to: %e \n", dlnns[s], ns[s]);
    }
    return; 
}

static void update_unknowns(double* S,double* dlnns,int nsp,double* ns,double* T, int verbose){
    /*
    Add corrections to unknown values (ignoring lagrange multipliers)
    Inputs:
        S : vector of corrections from matrix step [nel+1]
        dlnns : vector of species mole/mixture corrections [nsp]
        nsp : number of species
    Outputs:
        ns : vector of species mole/mixtures [nsp]
        T  : pointer to temperature guess (passed by reference!) [1]
    */
    int s;
    double lnns,lnT,n,lnn,lambda;
    const char pstring[] = "  s: %d lnns: % f rdlnns: % f dlnns: %f TR: % e lambda: % f\n"; 
    lnT = log(*T); // compute the log of the thing T is pointing to
    lambda = fmin(1.0, 0.5*fabs(lnT)/fabs(S[0]));
    *T = exp(lnT + lambda*S[0]); // thing pointed to by T set to exp(lnT + S[0]);

    n = 0.0; for (s=0; s<nsp; s++) n+=ns[s]; lnn=log(n);

    for (s=0; s<nsp; s++){
        if (ns[s]==0.0) {
            if (verbose>0) printf(pstring, s, -1.0/0.0, 0.0, dlnns[s], 0.0, 0.0);
            dlnns[s] = 0.0;
            continue;
        }
        lnns = log(ns[s]);
        lambda = fmin(1.0, fabs(lnn)/fabs(dlnns[s]));
        ns[s] = exp(lnns + lambda*dlnns[s]);

        if (ns[s]/n<TRACELIMIT){
            if (verbose>0) printf("    Locking species: %d (%f)\n", s, dlnns[s]);
            ns[s] = 0.0;
            dlnns[s] = 0.0; // This species is considered converged now
        }
    }
    return;
}

double temperature_guess(int nsp, double u, double M0, double* X0, double* lewis){
    /*
    Guess a first iteration temperature assuming constant Cv from 298 K
    Inputs:
        nsp   : Number of species
        u     : Target mixture internal energy (J/kg)
        M0    : Initial composition molecular weight (kg/mol)
        X0    : Intiial composition mole fractions [nsp]
        lewis : Nasa Lewis Thermodynamic Database Data [nsp*3*9]

    Output:
        T : Temperature Guess (K)
    */
    int s;
    double* lp;
    double uf,cv,T,ufs,cvs,Cps298,Hfs298,ns0;

    uf = 0.0;
    cv = 0.0;
    for (s=0; s<nsp; s++){
        lp = lewis + 9*3*s;
        Cps298 = compute_Cp0_R(298.15, lp)*Ru;
        Hfs298 = compute_H0_RT(298.15, lp)*Ru*298.15;

        ns0 = X0[s]/M0;
        ufs= ns0*(Cps298*298.15 - Hfs298);
        cvs= ns0*(Cps298 - Ru);

        uf += ufs;
        cv += cvs;
    }
    T = (u + uf)/cv;
    return T;
}

int solve_rhou(double rho,double u,double* X0,int nsp,int nel,double* lewis,double* M,double* a,
               double* X1, double* Teq, int verbose){
    /*
    Compute the equilibrium composition X1 at a fixed volume and internal energy 
    Inputs:
        rho   : target Density (kg/m3)
        u     : target internal energy (J/kg)
        X0    : Intial Mole fractions [nsp]
        nsp   : number of species 
        nel   : number of elements 
        lewis : Nasa Lewis Thermodynamic Database Data [nsp*3*9]
        M     : Molar Mass of each species (kg/mol) [nsp]
        a     : elemental composition array [nel,nsp]
        verbose: print debugging information

    Output:
        X1  : Equilibrium Mole Fraction [nsp]  
        Teq: Equilibrium Temperature 
    */
    double *A, *B, *S, *G0_RTs, *U0_RTs, *Cv0_Rs, *ns, *bi0, *dlnns, *errors; // Dynamic arrays
    int neq,s,i,k,errorcode,ntrace;
    double M0,n,M1,errorL2,errorL22,thing,T,errorrms;

    const double tol=1e-6;
    const int attempts=10;

    errorcode=0;
    neq= nel+1;
    A     = (double*) malloc(sizeof(double)*neq*neq); // Iteration Jacobian
    B     = (double*) malloc(sizeof(double)*neq);     // Iteration RHS
    S     = (double*) malloc(sizeof(double)*neq);     // Iteration unknown vector
    G0_RTs= (double*) malloc(sizeof(double)*nsp);     // Species Gibbs Free Energy
    U0_RTs= (double*) malloc(sizeof(double)*nsp);     // Species Internal Energy
    Cv0_Rs= (double*) malloc(sizeof(double)*nsp);     // Species Specific Heat @ Const. volume 
    ns    = (double*) malloc(sizeof(double)*nsp);     // Species moles/mixture mass
    bi0   = (double*) malloc(sizeof(double)*nel);     // starting composition coefficients
    dlnns = (double*) malloc(sizeof(double)*nsp);     // starting composition coefficients
    errors= (double*) malloc(sizeof(double)*nel);     // Nuclear constraint errors

    // Initialise Arrays and Iteration Guesses
    M0 = 0.0;
    for (s=0; s<nsp; s++) M0 += M[s]*X0[s];
    for (i=0; i<nel; i++){
        bi0[i] = 0.0;
        for (s=0; s<nsp; s++){
            bi0[i] += a[i*nsp + s]*X0[s]/M0;
            }
    }
    n = 0.0;
    for (s=0; s<nsp; s++) n += X0[s]/M0;
    for (s=0; s<nsp; s++) ns[s] = n/nsp;

    T = temperature_guess(nsp, u, M0, X0, lewis);
    if (verbose>0) printf("Guess T: %f\n", T);

    // Begin Iterations
    for (k=0; k<attempts; k++){
        Assemble_Matrices(a,bi0,rho,u,T,ns,nsp,nel,A,B,G0_RTs,U0_RTs,Cv0_Rs,lewis);
        errorcode = solve_matrix(A, B, S, neq);
        if (errorcode!=0) {
            handle_singularity(S,a,G0_RTs,U0_RTs,rho,T,ns,n,nsp,nel,dlnns,verbose);
            continue;
        }
        species_corrections(S,a,G0_RTs,U0_RTs,rho,T,ns,nsp,nel,dlnns,verbose);
        update_unknowns(S, dlnns, nsp, ns, &T, verbose);
        errorrms = constraint_errors(a, bi0, ns, nsp, nel, dlnns, verbose); //FIXME: Add T somehow?

        if (verbose>0){
            printf("iter %2d: [%f]",k,n);
            for (s=0; s<nsp; s++) printf(" %f",ns[s]);
            printf("  (%e)\n", errorrms);
        }

        // Exit loop if all but one species are trace 
        ntrace=0;
        for (s=0; s<nsp; s++) if (ns[s]==0.0) ntrace++;
        if (ntrace==nsp-1) { // Pseudo convergence criteria, all the species but one are trace
            for (s=0; s<nsp; s++) if (ns[s]!=0.0) i=s;
            n = 1.0/M[i];
            ns[i] = n;
            //errorL2 = 0.0;
            //errorL22= 0.0;
            errorrms = 0.0;
            if (verbose>0) printf("Pseudo convergence! Remaining species: (%d)\n", i);
        }

        // Exit loop if convergence is achieved
        if (errorrms<tol) break;

        // Exit loop if nans appearing in dlnns
        if (isnan(errorrms)) {
            printf("Solver nan'd, exiting!\n");
            errorcode=1;
            break;
        }

        if (k>=attempts) {
            printf("Solver not converged, exiting!\n");
            errorcode=1;
            break;
        }
    }
    
    if ((verbose>0)&&(errorcode==0)) printf("Converged in %d iter, error: %e\n", k, errorrms);
    // Compute output composition
    n = 0.0;
    for (s=0; s<nsp; s++) n += ns[s];
    M1 = 1.0/n;
    for (s=0; s<nsp; s++) X1[s] = M1*ns[s];
    *Teq = T;

    free(A);
    free(B);
    free(S);
    free(G0_RTs);
    free(U0_RTs);
    free(Cv0_Rs);
    free(ns);
    free(bi0);
    free(dlnns);
    free(errors);
    return errorcode;
}

#ifdef TEST
int main(){
    printf("Called main in rhou.c!\n");
    return 0;
}
#endif
