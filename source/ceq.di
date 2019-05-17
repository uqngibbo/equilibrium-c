module libceq;

extern(C) int pt(double p,double T,double* X0,int nsp,int nel,double* lewis,double* M,double* a,double* X1, int verbose);

extern(C) int rhou(double rho,double u,double* X0,int nsp,int nel,double* lewis,double* M,double* a,double* X1, double* T,int verbose);

extern(C) double get_u(double T, double* X, int nsp, double* lewis, double* M);

extern(C) int batch_pt(int N, double* p,double* T,double* X0,int nsp,int nel,double* lewis,double* M,double* a, double* X1, int verbose);

extern(C) int batch_rhou(int N, double* rho,double* u,double* X0,int nsp,int nel,double* lewis,double* M, double* a, double* X1, double* T, int verbose);

extern(C) int batch_u(int N, double* T, double* X, int nsp, double* lewis, double* M, double* u);
