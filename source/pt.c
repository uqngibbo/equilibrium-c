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
#include "pt.h"

static void constraint_errors(double* a, double* G0_RTs, double* bi0, double* pi, double p, double* ns, double n, int nsp, int nel, double* errors){
    /*
    Check the actual error in the Lagrange equations
    */ 
    int s,i;
    double lnn, lnp, bi;
    //lnn = log(n);
    //lnp = log(p/1e5); // Standard pressure for the tables is one BAR

    //for (s=0; s<nsp; s++){
    //    errors[s] = G0_RTs[s] + log(ns[s]) - lnn + lnp;
    //    for (i=0; i<nel; i++) errors[s] -= a[i*nsp + s]*pi[i];
    //}

    for (i=0; i<nel; i++) {
        bi = 0.0;
        for (s=0; s<nsp; s++) bi += a[i*nsp + s]*ns[s];
        errors[i] = bi - bi0[i];
    }
    return;
}

// static confines this function to this module
static void Assemble_Matrices(double* a,double* bi0,double* G0_RTs,double p,double* ns,double n,
                          int nsp,int nel,double* A, double* B){
    /*
    Construct Iteration Matrix for reduced Newton Rhapson step, (eqn 2.24 and 2.26 from cea_I)
    */
    double lnns, lnn, lnp, akjaijnj, akjnjmuj, mus_RTj, bk;
    double nss, nsmus, coeffsum;
    int k,neq,i,j,s;
    neq = nel+1;
    lnn = log(n);
    lnp = log(p/1e5); // Standard pressure for the tables is one BAR

    // Equation 2.24: k-> equation index, i-> variable index
    for (k=0; k<nel; k++){

        if (bi0[k]<1e-16) { // Check for missing missing element equations and Lock
            for (i=0; i<neq; i++) A[k*neq+i] = 0.0;
            A[k*neq + k+1] = 1.0;   
            B[k] = 0.0;
            continue;
        }

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
        B[k] = bi0[k] - bk + akjnjmuj;
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
    B[nel] = n - nss + nsmus;
    
    // Single Species Singularity Check: Set singular equations so that pii trivially equal to 1
    //for (i=0; i<neq; i++){
    //    coeffsum = 0.0; for (j=0; j<neq; j++) coeffsum += A[i*neq + j];
    //    if (coeffsum<1e-16){
    //        A[i*neq + i+1] = 1.0;   
    //        B[i] = 0.0;
    //    }
    //}

    for (i=0; i<neq; i++){
        printf("    [");
        for (j=0; j<neq; j++){
            printf("%f ", A[i*neq+j]);
        }
        printf("| %f]", B[i]);
        printf("\n");
    }
    return;
}

static void species_corrections(double* S,double* a,double* G0_RTs,double p,double n,double* ns,
                        int nsp, int nel, double* dlnns){
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
        dlnns[s] = -mu_RTs + dlnn + aispii;
        //printf("    dlnns[%d] = %f (%f %f %f)\n", s, dlnns[s], -mu_RTs, dlnn, aispii);
    }
    return; 
}

static void update_unknowns(double* S,double* dlnns,int nsp,int nel,double* ns,double* n, double* pi){
    /*
    Add corrections to unknown values (ignoring lagrange multipliers)
    Inputs:
        S : vector of corrections from matrix step [nel+1]
        dlnns : vector of species mole/mixture corrections [nsp]
        nsp : number of species
    Outputs:
        ns : vector of species mole/mixtures [nsp]
        n  : pointer to total moles/mixture (passed by reference!) [1]
        pi : Lagrange multipliers
    */
    int s,i;
    double lnns,lnn,n_copy,lambda;

    lnn = log(*n); // compute the log of the thing n is pointing to
    lambda = fmin(1.0, 0.5*fabs(lnn)/fabs(S[0]));
    n_copy = exp(lnn + lambda*S[0]); 
    *n = n_copy;   // thing pointed to by n set to exp(lnn + S[0]);
    for (i=0; i<nel; i++) pi[i] = S[i+1];

    for (s=0; s<nsp; s++){
        if (ns[s]==0.0) {
            //printf("    s: %d lnns: inf       dlnns: % f\n", s, dlnns[s]);
            continue;
        }
        lnns = log(ns[s]);

        lambda = fmin(1.0, 0.5*fabs(lnn)/fabs(dlnns[s]));
        //printf("    s: %d lnns: % f dlnns: % f lambda: % f\n", s, lnns, dlnns[s], lambda);
        ns[s] = exp(lnns + lambda*dlnns[s]);

        if (ns[s]/n_copy<TRACELIMIT){
            ns[s] = 0.0;
            dlnns[s] = 0.0; // This species is considered converged now
            //printf("    Locking species: %d\n", s);
        }
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
    double *A, *B, *S, *G0_RTs, *ns, *bi0, *dlnns, *errors, *pi; // Dynamic arrays
    double *lp;
    int neq,s,i,k,ntrace,errorcode,matrixerror;
    double M0,n,M1,errorL2,thing,errorL22;

    const double tol=1e-8;
    const int attempts=50;

    errorcode=0;
    matrixerror=0;
    neq= nel+1;
    A     = (double*) malloc(sizeof(double)*neq*neq); // Iteration Jacobian
    B     = (double*) malloc(sizeof(double)*neq);     // Iteration RHS
    S     = (double*) malloc(sizeof(double)*neq);     // Iteration unknown vector
    G0_RTs= (double*) malloc(sizeof(double)*nsp);     // Species Gibbs Free Energy
    ns    = (double*) malloc(sizeof(double)*nsp);     // Species moles/mixture mass
    bi0   = (double*) malloc(sizeof(double)*nel);     // starting composition coefficients
    dlnns = (double*) malloc(sizeof(double)*nsp);     // starting composition coefficients
    errors= (double*) malloc(sizeof(double)*nel);     // Trace Free energy constraint errors
    pi    = (double*) malloc(sizeof(double)*nel);     // Lagrange multipliers

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
    n*=1.1;

    for (s=0; s<nsp; s++){
        lp = lewis + 9*3*s;
        G0_RTs[s] = compute_G0_RT(T, lp);
    }

    // Auto lock species with missing elements
    for (i=0; i<nel; i++){
        if (bi0[i]<1e-16) {
            for (s=0; s<nsp; s++) {
                if (a[i*nsp + s]!=0) ns[s] = 0.0;
            }
        }
    }

    // Main Iteration Loop: 
    for (k=0; k<=attempts; k++){
        // 1: Perform an update of the equations
        Assemble_Matrices(a, bi0, G0_RTs, p, ns, n, nsp, nel, A, B);
        matrixerror = solve_matrix(A, B, S, neq);
        if (matrixerror!=0) {
             if (verbose) printf("    Singular Matrix!: Unlocking species\n");
             for (s=0; s<nsp; s++) ns[s] = fmax(2.0*n*TRACELIMIT, ns[s]); // Reset trace if singular
             continue;
        }
        species_corrections(S, a, G0_RTs, p, n, ns, nsp, nel, dlnns);
        update_unknowns(S, dlnns, nsp, nel, ns, &n, pi);
        constraint_errors(a, G0_RTs, bi0, pi, p, ns, n, nsp, nel, errors);

        // Compute remaining error by checking species corrections
        errorL2 = 0.0;
        for (s=0; s<nsp; s++) errorL2 += dlnns[s]*dlnns[s];
        errorL2 = pow(errorL2, 0.5);
        errorL22 = 0.0;
        for (i=0; i<nel; i++) errorL22 += errors[i]*errors[i];
        errorL22 = pow(errorL22, 0.5);

        if (verbose>0){
            printf("iter %2d: [%f]",k,n);
            for (s=0; s<nsp; s++) printf(" %f",ns[s]);
            printf("  (%e)\n", errorL2);
        }
        //if (verbose>0){
        //    printf("iter %d: [%e]",k,errorL2);
        //    for (i=0; i<nel; i++) printf(" %f",errors[i]);
        //    printf("  (%e)\n", errorL22);
        //}
        //if (verbose>0){
        //    printf("iter %2d: [%f]",k,n);
        //    for (i=0; i<nel; i++) printf(" %f",pi[i]);
        //    printf("  (%e)\n", errorL2);
        //}

        // Exit loop if all but one species are trace 
        ntrace=0;
        for (s=0; s<nsp; s++) if (ns[s]==0.0) ntrace++;
        if (ntrace==nsp-1) { // Pseudo convergence criteria, all the species but one are trace
            for (s=0; s<nsp; s++) if (ns[s]!=0.0) i=s;
            n = 1.0/M[i];
            ns[i] = n;
            errorL2 = 0.0;
            if (verbose>0) printf("Pseudo convergence! Remaining species: (%d)\n", i);
        }

        // Exit loop if convergence is achieved
        if (errorL2<tol) break;

        // Exit loop if nans appearing in dlnns
        if (isnan(errorL2)) {
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
    
    if ((verbose>0)&&(errorcode==0)) printf("Converged in %d iter, error: %e\n", k, errorL2);
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
    free(errors);
    free(pi);
    return errorcode;
}

#ifdef TEST
int main(){
    printf("Called main in pt.c!\n");
    return 0;
}
#endif
