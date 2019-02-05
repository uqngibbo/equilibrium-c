#ifndef ceq_h
#define ceq_h

int pt(double p,double T,double* X0,int nsp,int nel,double* lewis,double* M,double* a,double* X1,
       int verbose);

int rhou(double rho,double u,double* X0,int nsp,int nel,double* lewis,double* M,double* a,double* X1,
         double* T,int verbose);

int batch_pt(int N, double* p,double* T,double* X0,int nsp,int nel,double* lewis,double* M,double* a,
             double* X1, int verbose);

int batch_rhou(int N, double* rho,double* u,double* X0,int nsp,int nel,double* lewis,double* M,
               double* a, double* X1, double* T, int verbose);

#endif
