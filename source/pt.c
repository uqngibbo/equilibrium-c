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
#include "common.h"
#include "pt.h"

static void Assemble_Matrices(double* a,double* bi0,double* G0_RTs,double p,double* ns,double n,
                          int nsp,int nel,double* A, double* B){
    /*
    Construct Iteration Matrix for reduced Newton Rhapson step, (eqn 2.24 and 2.26 from cea_I)
    */
    double lnn, lnp, akjaijnj, akjnjmuj, mus_RTj, bk;
    double nss, nsmus;
    int k,neq,i,j,s;
    neq = nel+1;
    lnn = log(n);
    lnp = log(p/1e5); // Standard pressure for the tables is one BAR

    // Equation 2.24: k-> equation index, i-> variable index
    for (k=0; k<nel; k++){
        bk = 0.0; for (s=0; s<nsp; s++) bk += a[k*nsp + s]*ns[s];
        A[k*neq + 0] = bk;

        for (i=0; i<nel; i++){
            akjaijnj = 0.0;
            for (j=0; j<nsp; j++){
                akjaijnj += a[k*nsp+j]*a[i*nsp+j]*ns[j];
            }
            A[k*neq + i+1] = akjaijnj;
        }

        akjnjmuj = 0.0;
        for (j=0; j<nsp; j++){
            if (ns[j]==0.0) continue;
            mus_RTj = G0_RTs[j] + log(ns[j]) - lnn + lnp;
            akjnjmuj += a[k*nsp+j]*ns[j]*mus_RTj;

        }
        //B[k] = bi0[k] - bk + akjnjmuj;// CEA equations
        B[k] = bi0[k] + akjnjmuj;       // Newcode vers
        check_ill_posed_matrix_row(A, B, neq, k, 1);
    }

    // Equation 2.26 - > (only the pii entries, we're highjacking k here to go across the last row)
    for (k=0; k<nel; k++){
        bk = 0.0; for (s=0; s<nsp; s++) bk += a[k*nsp + s]*ns[s];
        A[nel*neq + k+1] = bk;
    }

    // Equation 2.26 - > (now the rest)
    nss = 0.0;
    nsmus = 0.0;
    for (j=0; j<nsp; j++){
        if (ns[j]==0.0) continue;
        mus_RTj = G0_RTs[j] + log(ns[j]) - lnn + lnp; // I guess its okay to compute this again
        nss += ns[j];
        nsmus += ns[j]*mus_RTj;
    }
    A[nel*neq + 0]  = nss - n;
    //B[nel] = n - nss + nsmus;  // CEA equations
    B[nel] = n + nsmus;          // Newcode vers
    
    //for (i=0; i<neq; i++){
    //    printf("    [");
    //    for (j=0; j<neq; j++){
    //        printf("%f, ", A[i*neq+j]);
    //    }
    //    printf("  %f]", B[i]);
    //    printf("\n");
    //}
    return;
}

static double compute_residual(double* pij, double* a, double* G0_RTs, double p, double T, double n,
                               double* ns, double* bi0, int nsp, int nel, int verbose){
    /*
    Compute the L2 of the PT Lagrangian derivatives. Note as per the paper, we actually compute
    ns[s] times equation (44). This ensures that the equations are nonsingular for ns[s] -> 0.0
    but still have the same root as the original equations.

    Inputs:
        pij    : Corrections (dlog(n), pi1, pi2, pi3 ...)  [nel+1]
        a      : elemental composition array [nel,nsp]
        G0_RTs : Gibbs free energy of each species, divided by RT [nsp]
        p      : pressure  (Pa)
        T      : Temperature (K)
        n      : total moles/mixture kg
        ns     : species moles/mixture kg [nsp]
        nsp    : total number of species
        nel    : total number of elements

    Outputs:
        double : L2 norm of the residual of the nonlinear equations being solved.
    */

    double residual = 0.0;

    for (int s=0; s<nsp; s++){
        // ns*log(ns) is badly behaved at ns==0.0, but should be fine at any other nonnegative number
        double nss = ns[s];
        double nslogns;
        if (nss==0.0) { nslogns=0.0; } else { nslogns = nss*log(nss); }

        // Equation ?? from the eqc paper
        //double Rs = nss*G0_RTs[s] + nslogns - nss*log(n) + nss*log(p/1e5);         // CEA equations
        double Rs = nss*G0_RTs[s] + nslogns - nss*log(n) + nss*log(p/1e5) + 1.0*nss; // Newcode vers
        double pijj = 0.0;
        for (int j=0; j<nel; j++){
            pijj += pij[1+j]*a[j*nsp + s];
        }
        Rs -= nss*pijj;
        if (verbose>1) printf("ns*Rss[%d]=%e\n", s, Rs);
        residual += Rs*Rs;
    }

    // Equation ?? from the eqc paper
    for (int j=0; j<nel; j++){
        double Rs = 0.0;
        for (int s=0; s<nsp; s++){
            Rs += a[j*nsp + s]*ns[s];
        }
        Rs -= bi0[j];
        if (verbose>1) printf("Rse[%d]=%e\n", j, Rs);
        residual += Rs*Rs;
    }

    // Equation ?? from the eqc paper
    double Rs = n;
    for (int s=0; s<nsp; s++){
        Rs -= ns[s];
    }
    if (verbose>1) printf("Rsn=%e\n", Rs);
    residual += Rs*Rs;

    return sqrt(residual);
}

static double compute_lagrangian(double* pij, double* a, double* G0_RTs, double p, double T, double n,
                                 double* ns, double* bi0, int nsp, int nel, int verbose){

    double L = 0.0;

    for (int s=0; s<nsp; s++){
        // ns*log(ns) is badly behaved at ns==0.0, but should be fine at any other nonnegative number
        double nss = ns[s];
        double nslogns;
        if (nss==0.0) { nslogns=0.0; } else { nslogns = nss*log(nss); }

        L += Ru*T*(nss*G0_RTs[s] + nslogns - nss*log(n) + nss*log(p/1e5));
    }

    for (int j=0; j<nel; j++){
        double ajsns = 0.0; for (int s=0; s<nsp; s++) ajsns += a[j*nsp + s]*ns[s];
        double lambda_j = -1.0*Ru*T*pij[j+1];
        L += lambda_j*(ajsns - bi0[j]);
    }

    return L;
}

static void species_corrections(double* S,double* a,double* G0_RTs,double p,double n,double* ns,
                        int nsp, int nel, double* dlnns, int verbose){
    /*
    Compute delta_log(ns) from the reduced iteration equations from 
    equation 2.18m using the other deltas in S
    Inputs:
        S      : Corrections (pi1, pi2, pi3 ... dlog(n) [nel+1]
        a      : elemental composition array [nel,nsp]
        G0_RTs : Gibbs free energy of each species, divided by RT [nsp]
        p      : pressure 
        n      : total moles/mixture kg 
        ns     : species moles/mixture kg [nsp]
        nsp    : total number of species
        nel    : total  number of elements 

    Outputs:
        dllns : change in log(ns) [nsp]
    */
    double dlnn,aispii,mu_RTs,lnn,lnp;
    int s,i;
    dlnn = S[0];
    lnn = log(n);
    lnp = log(p/1e5);

    for (s=0; s<nsp; s++) {
        if (ns[s]==0.0) { dlnns[s] = 0.0; continue;}
        mu_RTs = G0_RTs[s] + log(ns[s]) - lnn + lnp;

        aispii = 0.0;
        for (i=0; i<nel; i++){
            aispii += a[i*nsp+s]*S[i+1]; // S[i+1] = pi_i, the lagrange multiplier
        }
        //dlnns[s] = -mu_RTs + dlnn + aispii;     // CEA equations
        dlnns[s] = -mu_RTs + dlnn + aispii - 1.0; // Newcode vers
        //printf("    dlnns[%d] = %f (%f %f %f)\n", s, dlnns[s], -mu_RTs, dlnn, aispii);
    }
    return; 
}

static void handle_singularity(double* S,double* a,double* G0_RTs,double p,double n,double* ns,
                        int nsp, int nel, double* dlnns, int verbose){
    /*
    Handle crash by possibly resetting species compositions to fix
    Inputs:
        S      : Corrections array (pi1, pi2, pi3 ... dlog(n) [nel+1]
        a      : elemental composition array [nel,nsp]
        G0_RTs : Gibbs free energy of each species, divided by RT [nsp]
        p      : pressure 
        n      : total moles/mixture kg 
        ns     : species moles/mixture kg [nsp]
        nsp    : total number of species
        nel    : total  number of elements 
        verbose: flag to print debugging  information
    */
    const double RESET=4.0;
    int s;

    if (verbose>0) printf("    ### Singular Matrix!: Unlocking to %f\n",log(RESET*n*TRACELIMIT));

    for (s=0; s<nsp; s++){
        if (ns[s]!=0.0) continue;  // Ignore non trace species
        ns[s] = fmax(RESET*n*TRACELIMIT, ns[s]); // Reset trace trace species to a small but finite number
        species_corrections(S, a, G0_RTs, p, n, ns, nsp, nel, dlnns, 0); // approximately predict dlnns
        if (dlnns[s]<0.0) ns[s] = 0.0; // Re-zero any species with negative predicted dlnns
        if (verbose>1) printf("   faux dlnns: %f changed to: %e \n", dlnns[s], ns[s]);
    }
    return; 
}

static void update_unknowns(double* S,double* dlnns,int nsp,int nel,double* ns,double* n,int verbose){
    /*
    Add corrections to unknown values (ignoring lagrange multipliers)
    Inputs:
        S : vector of corrections from matrix step [nel+1]
        dlnns : vector of species mole/mixture corrections [nsp]
        nsp : number of species
    Outputs:
        ns : vector of species mole/mixtures [nsp]
        n  : pointer to total moles/mixture (passed by reference!) [1]
    */
    int s;
    double lnns,lnn,n_copy,lambda,newns,rdlnns;
    const char pstring[] = "  s: %d lnns: % f rdlnns: % f dlnns: %f TR: % e lambda: % f\n"; 

    lnn = log(*n); // compute the log of the thing n is pointing to
    lambda = update_limit_factor(lnn, S[0], relaxation_limit);
    n_copy = exp(lnn + lambda*S[0]); 
    *n = n_copy;   // thing pointed to by n set to exp(lnn + S[0]);

    for (s=0; s<nsp; s++){
        if (ns[s]==0.0) {
            if (verbose>1) printf(pstring, s, 0.0, 0.0, dlnns[s], 0.0, 0.0);
            dlnns[s] = 0.0;
            continue;
        }
        lnns = log(ns[s]);
        lambda = update_limit_factor(lnn, dlnns[s], relaxation_limit);
        newns = exp(lnns + lambda*dlnns[s]);
        rdlnns = log(newns) - lnns;
        ns[s] = newns;
        lnns = log(newns);

        if (verbose>1) printf(pstring, s, lnns, rdlnns, dlnns[s], ns[s]/n_copy/TRACELIMIT, lambda);
    }

    return;
}

int solve_pt(double p,double T,double* X0,int nsp,int nel,double* lewis,double* M,double* a,double* X1,int verbose){
    /*
    Compute the equilibrium composition X1 at a fixed temperature and pressure
    Inputs:
        p     : Pressure (Pa)
        T     : Temperature (K)
        X0    : Intial Mole fractions [nsp]
        nsp   : number of species 
        nel   : number of elements 
        lewis : Nasa Lewis Thermodynamic Database Data [nsp*3*9]
        M     : Molar Mass of each species (kg/mol) [nsp]
        a     : elemental composition array [nel,nsp]
        verbose: print debugging information

    Output:
        X1 : Equilibrium Mole Fraction [nsp]  
    */
    double *A, *B, *S, *G0_RTs, *ns, *bi0, *dlnns; // Dynamic arrays
    double *lp;
    int neq,s,i,k,errorcode;
    double n,M1,errorrms;

    errorcode=0;
    neq= nel+1;
    errorrms=1e99;
    A     = (double*) malloc(sizeof(double)*neq*neq); // Iteration Jacobian
    B     = (double*) malloc(sizeof(double)*neq);     // Iteration RHS
    S     = (double*) malloc(sizeof(double)*neq);     // Iteration unknown vector
    G0_RTs= (double*) malloc(sizeof(double)*nsp);     // Species Gibbs Free Energy
    ns    = (double*) malloc(sizeof(double)*nsp);     // Species moles/mixture mass
    bi0   = (double*) malloc(sizeof(double)*nel);     // starting composition coefficients
    dlnns = (double*) malloc(sizeof(double)*nsp);     // starting composition coefficients

    composition_guess(a, M, X0, nsp, nel, ns, &n, bi0);
    n*=1.1;

    for (s=0; s<nsp; s++){
        lp = lewis + 9*3*s;
        G0_RTs[s] = compute_G0_RT(T, lp);
    }

    // Main Iteration Loop: 
    for (k=0; k<=attempts; k++){
        // 1: Perform an update of the equations
        Assemble_Matrices(a, bi0, G0_RTs, p, ns, n, nsp, nel, A, B);
        errorcode = solve_matrix(A, B, S, neq);
        if (errorcode!=0) {
            handle_singularity(S, a, G0_RTs, p, n, ns, nsp, nel, dlnns, verbose);
            continue;
        }
        species_corrections(S, a, G0_RTs, p, n, ns, nsp, nel, dlnns, verbose);
        update_unknowns(S, dlnns, nsp, nel, ns, &n, verbose);
        handle_trace_species_locking(a, n, nsp, nel, ns, bi0, dlnns, verbose);
        errorrms = compute_residual(S, a, G0_RTs, p, T, n, ns, bi0, nsp, nel, verbose);


        if (verbose>0){
            printf("iter %2d: [%f]",k,n);
            for (s=0; s<nsp; s++) printf(" %f",ns[s]);
            printf("  (%e)\n", errorrms);
        }


        // Exit loop if all but one species are trace 
        i = all_but_one_species_are_trace(nsp, ns);
        if (i>=0){
            if (verbose>0) printf("Pseudo convergence! Remaining species: (%d)\n", i);
            n = 1.0/M[i];
            ns[i] = n;
            errorrms = 0.0;
            break;
        }

        // Exit loop if convergence is achieved
        if (errorrms<tol) break;

        // Exit loop if nans appearing in dlnns
        if (isnan(errorrms)) {
            printf("Solver nan'd, exiting!\n");
            errorcode=1;
            break;
        }

        // Exit loop if too many attempts are undertaken
        if (k==attempts) {
            printf("Solver not converged, exiting!\n");
            errorcode=1;
            break;
        }
    }
    
    if ((verbose>0)&&(errorcode==0)) printf("Converged in %d iter, error: %e\n", k, errorrms);
    if ((verbose>0)&&(errorcode>0)) printf("Convergence failure in %d iter, error: %e\n", k, errorrms);

    // For verification purposes, we might want to numerically differentiate the
    // Lagrangian to check that we have correctly found the actual stationary point.
    if (verbose>0){
        double eps = 1e-7;
        double L = compute_lagrangian(S, a, G0_RTs, p, T, n, ns, bi0, nsp, nel, verbose);
        printf("Lagrangian:[%e]\n   ",L);
        for (s=0; s<nsp; s++) {
            double ns_save = ns[s];
            double perturb = ns[s]*eps;
            ns[s] += perturb;
            double L2 = compute_lagrangian(S, a, G0_RTs, p, T, n, ns, bi0, nsp, nel, verbose);
            ns[s] = ns_save;
            double dLdns = (L2-L)/perturb;
            printf(" dLdn[%d]= %e",s, dLdns);
        }
        printf("\n");
    }
    // Compute output composition
    M1 = 1.0/n;
    for (s=0; s<nsp; s++) X1[s] = M1*ns[s];

    free(A);
    free(B);
    free(S);
    free(G0_RTs);
    free(ns);
    free(bi0);
    free(dlnns);
    return errorcode;
}

#ifdef TEST
int main(){
    printf("Called main in pt.c!\n");
    return 0;
}
#endif
