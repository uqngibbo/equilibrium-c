#ifndef ceq_h
#define ceq_h

void pt_Assemble_Matrices(double* a,double* bi0,double* G0_RTs,double p,double* ns,double n,
                          int nsp,int nel,double* A, double* B);
void species_corrections(double* S,double* a,double* G0_RTs,double p,double n,double* ns,
                        int nel, int nsp, double* dlnns);
void update_unknowns(double* S,double* dlnns,int nsp,double* ns,double* n);
int pt(double p,double T,double* X0,int nsp,int nel,double* lewis,double* M,double* a,double* X1);

#endif
