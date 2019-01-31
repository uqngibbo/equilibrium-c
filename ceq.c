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

#include "pt.h"
#include "rhou.h"
#include "ceq.h"


int pt(double p,double T,double* X0,int nsp,int nel,double* lewis,double* M,double* a,double* X1){
    /*
    Compute the equilibrium composition X1 at a fixed temperature and pressure (Wrapper Function)
    */
    return solve_pt(p, T, X0, nsp, nel, lewis, M, a, X1);
}

int rhou(double rho,double u,double* X0,int nsp,int nel,double* lewis,double* M,double* a,
         double* X1){
    /*
    Compute the equilibrium composition X1 at a fixed density and internal energy (Wrapper Function)
    */
    return solve_rhou(rho, u, X0,nsp,nel, lewis, M, a, X1);
}


int main(){
    printf("Called main in ceq.c!\n");
    return 0;
}
