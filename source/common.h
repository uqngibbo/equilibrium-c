#ifndef common_h
#define common_h

int all_but_one_species_are_trace(int nsp, double* ns);
double constraint_errors(double* S, double* a,double* bi0,double* ns,int nsp,int nel,double* dlnns,int verbose);
void composition_guess(double* a,double* M,double* X0,int nsp,int nel,double* ns,double* np,double* bi0);

#endif
