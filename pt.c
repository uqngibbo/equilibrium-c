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

// static confines this function to this module
static void Assemble_Matrices(double* a,double* bi0,double* G0_RTs,double p,double* ns,double n,
                          int nsp,int nel,double* A, double* B){
    /*
    Construct Iteration Matrix for reduced Newton Rhapson step, (eqn 2.24 and 2.26 from cea_I)
    */
    double lnns, lnn, lnp, akjaijnj, akjnjmuj, mus_RTj, bk;
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
            mus_RTj = G0_RTs[j] + log(ns[j]) - lnn + lnp;
            akjnjmuj += a[k*nsp+j]*ns[j]*mus_RTj;

        }
        B[k] = bi0[k] - bk + akjnjmuj;

        // Equation 2.26 - > (only the pii entries, we're highjacking k here to go across the last row)
        A[nel*neq + k+1] = bk;
    }
    // Equation 2.26 - > (now the rest)
    nss = 0.0;
    nsmus = 0.0;
    for (j=0; j<nsp; j++){
        mus_RTj = G0_RTs[j] + log(ns[j]) - lnn + lnp; // I guess its okay to compute this again
        nss += ns[j];
        nsmus += ns[j]*mus_RTj;
    }
    A[nel*neq + 0]  = nss - n;
    B[nel] = n - nss + nsmus;

    //for (i=0; i<neq; i++){
    //    for (j=0; j<neq; j++){
    //        printf("%f ", A[i*neq+j]);
    //    }
    //    printf("| %f\n", B[i]);
    //}
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
        mu_RTs = G0_RTs[s] + log(ns[s]) - lnn + lnp;

        aispii = 0.0;
        for (i=0; i<nel; i++){
            aispii += a[i*nsp+s]*S[i+1]; // S[i+1] = pi_i, the lagrange multiplier
        }
        dlnns[s] = -mu_RTs + dlnn + aispii;
    }
    return; 
}

static void update_unknowns(double* S,double* dlnns,int nsp,double* ns,double* n){
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
    double lnns,lnn;
    for (s=0; s<nsp; s++){
        lnns = log(ns[s]);
        ns[s] = exp(lnns + dlnns[s]);
    }
    lnn = log(*n); // compute the log of the thing n is pointing to
    *n = exp(lnn + S[0]); // thing pointed to by n set to exp(lnn + S[0]);
    return;
}

int solve_pt(double p,double T,double* X0,int nsp,int nel,double* lewis,double* M,double* a,double* X1){
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

    Output:
        X1 : Equilibrium Mole Fraction [nsp]  
    */
    double *A, *B, *S, *G0_RTs, *ns, *bi0, *dlnns; // Dynamic arrays
    double *lp;
    int neq,s,i,k;
    double M0,n,M1,errorL2,thing;

    const double tol=1e-6;
    const int attempts=10;

    neq= nel+1;
    A     = (double*) malloc(sizeof(double)*neq*neq); // Iteration Jacobian
    B     = (double*) malloc(sizeof(double)*neq);     // Iteration RHS
    S     = (double*) malloc(sizeof(double)*neq);     // Iteration unknown vector
    G0_RTs= (double*) malloc(sizeof(double)*nsp);     // Species Gibbs Free Energy
    ns    = (double*) malloc(sizeof(double)*nsp);     // Species moles/mixture mass
    bi0   = (double*) malloc(sizeof(double)*nel);     // starting composition coefficients
    dlnns = (double*) malloc(sizeof(double)*nsp);     // starting composition coefficients

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

    // Begin Iterations
    for (k=0; k<attempts; k++){
        Assemble_Matrices(a, bi0, G0_RTs, p, ns, n, nsp, nel, A, B);
        solve_matrix(A, B, S, neq);
        species_corrections(S, a, G0_RTs, p, n, ns, nsp, nel, dlnns);
        update_unknowns(S, dlnns, nsp, ns, &n);

        errorL2 = 0.0;
        for (s=0; s<nsp; s++) errorL2 += dlnns[s]*dlnns[s];
        errorL2 = pow(errorL2, 0.5);
        if (errorL2<tol) break;

        if (k>=attempts) {
            printf("Solver not converged, exiting!\n");
            return 1;
        }
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
    return 0;
}

#ifdef TEST
int main(){
    printf("Called main in pt.c!\n");
    return 0;
}
#endif
