/*
C library for equilibrium chemistry calculations: linalg module

@author: Nick Gibbons
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "linalg.h"

void solve_matrix(double* A, double* B, double *X, int N){
    /*
    Bare bones Gaussian matrix solution A*X = B 
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
    //for (i=0; i<N; i++) {
    //    printf("[");
    //    for (j=0; j<N; j++){
    //        printf("%f ", A[i*N+j]);
    //    }
    //    printf("]\n");
    //}
    //printf("\n");

    } // end while loop for row reduction

    // Now compute X using back substitution
    for (k=N-1; k>=0; k--){
        a = 0.0;
        for (j=k+1; j<N; j++) a+= A[k*N+j]*X[j]; // Should skip when k=N
        X[k] = B[k] - a;
    }
    return;
}

void test_solve_matrix(){
    double XX[4];
    int N=4;
    int i;

    double A[16] = {0, 3, 5, 4, 6, 1, 7, 5, 1, 3, 2, 6, 6, 0, 2, 0};
    double X[4] = {6.0, 1.0, 4.0, 1.0};
    double B[4] = {27.0, 70.0, 23.0, 44.0};

    solve_matrix(A, B, XX, N);

    for (i=0; i<N; i++){
        printf("%d X: %f XX: %f\n",i , X[i], XX[i]);
    }
    return;
}

int main(){
    test_solve_matrix();
    return 0;
}
