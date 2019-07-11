/*
C library for equilibrium chemistry calculations
 - This file contains useful utility functions shared between pt.c and rhou.c

@author: Nick Gibbons
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void composition_guess(double* a,double* M,double* X0,int nsp,int nel,double* ns,double* np,double* bi0){
    /*
    Unified setup of initial composition from mole fractions
    Inputs:
        a      : elemental composition array [nel,nsp]
        X0    : Intial Mole fractions [nsp]
        M      : species molecular masses [nsp]
        nsp    : total number of species
        nel    : total  number of elements 
    Outputs:
        ns     : species moles/mixture kg [nsp]
        n      : total moles/mixture kg  [1]
        bi0    : Initial Nuclear moles/mixture [nel]
    */ 
    int i,s;
    double M0,n;

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
    *np = n;

    // Auto lock species with missing elements
    for (i=0; i<nel; i++){
        if (bi0[i]<1e-16) {
            for (s=0; s<nsp; s++) {
                if (a[i*nsp + s]!=0) ns[s] = 0.0;
            }
        }
    }
    return;
}

int all_but_one_species_are_trace(int nsp, double* ns){
    /*
    If only one species is left in the calculation assume we have found the answer
    Inputs:
        nsp : number of species 
        ns  : species moles/mixture kg [nsp]
    Returns:
        -1 if false, equal to the index of the nontrace species otherwise
    */
    int i,s,ntrace;

    ntrace=0;
    for (s=0; s<nsp; s++) if (ns[s]==0.0) ntrace++;

    if (ntrace==nsp-1) { // Pseudo convergence criteria, all the species but one are trace
        for (s=0; s<nsp; s++) if (ns[s]!=0.0) i=s;
        return i;
    }
    else {
        return -1;
    }
}

double constraint_errors(double* S, double* a,double* bi0,double* ns,int nsp,int nel,double* dlnns,int verbose){
    /*
    Unified computation of current error to determine whether to break loop
    Inputs:
        S      : Corrections  [nel+1]
        a      : elemental composition array [nel,nsp]
        bi0    : Initial Nuclear moles/mixture [nel]
        ns     : species moles/mixture kg [nsp]
        nsp    : total number of species
        nel    : total  number of elements 
        dllns  : raw change in log(ns) [nsp]
        verbose: verbose flag to print debugging information

    Outputs:
        unified error number (method as determined by this function)
    */ 
    int s,i;
    double bi,errorrms,error;

    // Compute the change in current variables (note this is the unlimited dlnns! not the real change)
    errorrms = 0.0;

    errorrms += S[0]*S[0]; // Either dlnT or dlnn

    for (s=0; s<nsp; s++) errorrms += dlnns[s]*dlnns[s];

    for (i=0; i<nel; i++) {
        bi = 0.0;
        for (s=0; s<nsp; s++) bi += a[i*nsp + s]*ns[s];
        error = bi - bi0[i];
        errorrms += error*error;
    }
    errorrms /= (nsp+nel+1);
    errorrms = sqrt(errorrms);
    return errorrms;
}
