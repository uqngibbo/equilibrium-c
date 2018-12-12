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
#include "ceq.h"

void pt_Assemble_Matrices(double* a,double* bi0,double* G0_RTs,double p,double* ns,double n,
                          int nsp,int nel,double* A, double* B){
    /*
    Construct Iteration Matrix for reduced Newton Rhapson step, (eqn 2.24 and 2.26 from cea_I)
    */
    double lnns, lnn, lnp, akjaijnj, akjnjmuj, mus, nsmus;
    int k,neq,i,j;
    neq = nel+1;
    lnn = log(n);
    lnp = log(p/1e5);

    for (k=0; k<nel; k++){
        bk = 0.0; for (s=0; s<nsp; s++) bk += a[k*nel + s]*ns[s];
        A[k*neq + 0] = bk;

        for (i=0; i<nel; i++){
            akjaijnj = 0.0;
            for (j=0; j<nsp; j++){
                akjaijnj += a[k*nel+j]*a[i*nel+j]*ns[j];
            }
            A[k*neq + i+1] = akjaingn;
        }
        akjnjmuj = 0.0;
        for (j=0; j<nsp; j++){

        }


    }

}

int pt(double p,double T,double* X0,int nsp,int nel,double* lewis,double* M,double* a,double* X1){
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
    double* A, B, S, G0_RTs, pii, ns, bi0, dlnns; // Dynamic arrays
    double* lp;
    int neq,s,i,k;
    double M0,n,M1;

    const double tol=1e-6;

    neq= nel+1;
    A     = (double*) malloc(sizeof(double)*neq*neq); // Iteration Jacobian
    B     = (double*) malloc(sizeof(double)*neq);     // Iteration RHS
    S     = (double*) malloc(sizeof(double)*neq);     // Iteration unknown vector
    G0_RTs= (double*) malloc(sizeof(double)*nsp);     // Species Gibbs Free Energy
    pii   = (double*) malloc(sizeof(double)*nel);     // Lagrange Multipliers
    ns    = (double*) malloc(sizeof(double)*nsp);     // Species moles/mixture mass
    bi0   = (double*) malloc(sizeof(double)*nel);     // starting composition coefficients
    dlnns = (double*) malloc(sizeof(double)*nsp);     // starting composition coefficients

    // Initialise Arrays and Iteration Guesses
    M0 = 0.0;
    for (s=0; s<nsp; s++) M0 += M[s]*X0[s];
    for (i=0; i<nel; i++){
        bi0[i] = 0.0;
        for (s=0; s<nsp; s++) bi0[i] += a[i*nel + s]*X0[s]/M0;
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
    for (k=0; k<20; k++){
        pt_Assemble_Matrices(a, bi0, G0_RTs, p, ns, n, A, B);
        solve_matrix(A, B, S, neq);
        species_corrections(S, a, bi0, G0_RTs, p, n, ns, dlnns);
        update_unknowns(S, dlnns, ns, &n, pi);

        errorL2 = 0.0;
        for (s=0; s<nsp; s++) errorL2 += dlnns[s]*dlnns[s];
        errorL2 = pow(errorL2, 0.5);
        if (errorL2<tol) break;

        if (k>19) {
            printf("Solver not converged, exiting!");
            return 1;
        }
    }
    
    // Compute output composition
    M1 = 1.0/n;
    for (s=0; s<nsp; s++) X1[s] = M1*ns[s];

    free(A);
    free(B);
    free(S);
    free(G0_RT);
    free(pii);
    free(ns);
    free(bi0);
    free(dlnns);
    return 0;
}

void test_pt(){

}

int main(){
    return 0;
}
