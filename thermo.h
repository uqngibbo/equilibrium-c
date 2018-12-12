#ifndef thermo_h
#define thermo_h

const double Ru=8.3144621;      // Universal Gas Constant (J/mol/K)
double compute_H0_RT(double Tin, double* lewis);
double compute_S0_R(double Tin, double* lewis);
double compute_G0_RT(double T, double* lewis);

#endif
