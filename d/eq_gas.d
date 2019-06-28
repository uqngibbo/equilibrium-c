/**
 * eq_gas.d
 * Equilibrium gas mixture based on ceq library
 *
 * Author: Nick Gibbons 
 * Notes: - Solving works
 *        - Consider issues with initial concentrations
 *        - Consider issues with small species getting "deadlocked" into Ys=0
 *        - Consider OOP reimplementation if more solver modes are needed (rhoT etc)
 */

module gas.eq_gas;

import std.math;
import std.stdio;
import std.string;
import std.conv : to;
import util.lua;
import util.lua_service;
import std.algorithm: canFind, countUntil;
import std.algorithm.sorting: sort;
import std.range: enumerate, take;
import nm.number;
import nm.bbla;
import gas.thermo.cea_thermo_curves;
import gas.gas_model;
import gas.gas_state;
import gas.therm_perf_gas;

//FIXME: Remove hardcoded DBPATH and import physical constants from the d module
const string DBPATH = "/home/qungibbo/programs/ceq/thermo.inp";
const double Ru=8.3144621;      // Universal Gas Constant (J/mol/K)
const double TRACELIMIT=1e-3;   // Trace species limiter (for ns/n)
const double P_atm=101.325e3;

class EqGas: ThermallyPerfectGas {
public:
    this(lua_State* L) {
        super(L);

        // We're going to trick the flow solver into thinking that we have only one species
        nsp = to!int(_n_species);
        species_list.length = nsp;
        M.length = nsp;
        species_list = _species_names.dup();
        M = _mol_masses.dup();

        double[string] species_elements;
        double[string][] element_map;
        string[] element_set;

        foreach(species; species_list){
            species_elements.clear();
            lua_getglobal(L, "db");
            lua_getfield(L, -1, species.toStringz);
            lua_getfield(L, -1, "atomicConstituents");

            // Why does this completely idiomatic operation take TEN lines ...
            // Why do I have to micromanage the key numbers in weird inconsistent ways
            // Why does lua next automatically mess with the stack....
            // WHYYYYYYYYYYYYY
            lua_pushnil(L); // dummy first key
            while (lua_next(L, -2) != 0) { // -1 is the dummy key, -2 is the atomicConstituents table
                // lua_next removes the dummy key and puts the first key value pair on the stack
                // the key at -2, value at -1
                string key = to!string(lua_tostring(L, -2));
                int value = to!int(lua_tointeger(L, -1));
                species_elements[key] = to!double(value);
                lua_pop(L, 1); // discard value but keep key so that lua_next can remove it (?!)
            }
            lua_pop(L, 1); // remove atomicConstituents (lua_next removed the key when it broke loop)
            lua_pop(L, 1); // remove species table
            lua_pop(L, 1); // remove db table
            element_map ~= species_elements.dup();

        }

        compile_element_set(element_map, element_set);
        compile_element_matrix(species_list, element_map, element_set, a);

        nel = to!int(element_set.length);
        X0.length = nsp; X1.length = nsp;

        neq= nel+1;
        S.length      = neq;     // Iteration unknown vector
        G0_RTs.length = nsp;     // Species Gibbs Free Energy
        U0_RTs.length = nsp;     // Species Internal Energys per mole
        Cv0_Rs.length = nsp;     // Species Molar specific heat at constant volume
        ns.length     = nsp;     // Species moles/mixture mass
        bi0.length    = nel;     // starting composition coefficients
        dlnns.length  = nsp;     // starting composition coefficients
        AB = new Matrix!double(neq, neq+1);
        
    } // End Lua state contructor

    this(in string fname)
    {
        auto L = init_lua_State();
        doLuaFile(L, fname);
        this(L);
        lua_close(L); // We no longer need the Lua interpreter.
    } // end constructor from a Lua file

    override string toString() const
    {
        return "EqGas(species=[TODO])";
    }

    override void update_thermo_from_pT(GasState Q) 
    {
        int error;
        //printf("Called update_thermo_from_pt\n");
        //printf("Q.T: %f  Q.p: %f Q.rho: %f Q.u: %f\n", Q.T, Q.p, Q.rho, Q.u);
        //printf("Q.YCO2: %f  Q.YCO: %f Q.YO2: %f\n", Q.massf[0], Q.massf[1], Q.massf[2]);

        massf2molef(Q, X0); 
        error = solve_pt(Q.p,Q.T,X0,M,a,X1,0);
        molef2massf(X1, Q);
        _pgMixEOS.update_density(Q);
        _tpgMixEOS.update_energy(Q);

    }
    override void update_thermo_from_rhou(GasState Q)
    {
        int error;
        printf("Called update_thermo_from_rhou\n");
        printf("Q.T: %f  Q.p: %f Q.rho: %f Q.u: %f\n", Q.T, Q.p, Q.rho, Q.u);
        printf("Q.YCO2: %f  Q.YCO: %f Q.YO2: %f\n", Q.massf[0], Q.massf[1], Q.massf[2]);

        massf2molef(Q, X0); 
        error = solve_rhou(Q.rho,Q.u,X0,M,a,X1,Q.T,1);
        molef2massf(X1, Q);
        _pgMixEOS.update_pressure(Q);
    }

    override void update_thermo_from_rhoT(GasState Q) const
    {
        throw new Exception("EqGas update_thermo_from_rhoT not implemented.");
    }
    override void update_thermo_from_rhop(GasState Q) const
    {
        throw new Exception("EqGas update_thermo_from_rhop not implemented.");
    }
    
    override void update_thermo_from_ps(GasState Q, number s) const
    {
        throw new Exception("EqGas update_thermo_from_ps not implemented.");
    }
    override void update_thermo_from_hs(GasState Q, number h, number s) const
    {
        throw new Exception("EqGas update_thermo_from_hs not implemented.");
    }

    // In case another part of the code tries to ask for a specific species, give them the mass average
    override number enthalpy(in GasState Q, int isp)
    {
        int s;
        for (s=0; s<nsp; s++) {
            _h[s] = _curves[s].eval_h(Q.T);
        }
        return mass_average(Q, _h);
    }
    override number entropy(in GasState Q, int isp)
    {
        int s;
        for (s=0; s<nsp; s++) {
            _s[s] = _curves[s].eval_s(Q.T) - _R[s]*log(Q.p/1e5);
        }
        return mass_average(Q, _s);
    }

    // FIXME: Possible issue with species names property 
    @nogc @property override uint n_Conserved_species() const { return 0; } // FIXME: possibly make it 1?

private:
    int nel, nsp, neq;
    string[] species_list;
    double[] a;
    double[] M;
    double[] X0,X1;
    double[] S, G0_RTs, U0_RTs, Cv0_Rs, ns, bi0, dlnns; 
    Matrix!double AB;

    void compile_element_set(double[string][] elements_array, ref string[] elements){
        /*
        Get an alphabetical list of all the elements in the complete set of species
        Inputs: 
             elements_array : array of associative arrays mapping species to elements
        Outputs:
            elements : list of elements in the entire system
        */
    
        foreach(e; elements_array){
            foreach(key; e.keys()){
                if (!elements.canFind(key)) elements ~= key;
            }
        }
        elements.sort();
        return;
    }
    
    void compile_element_matrix(string[] speciesList, double[string][] elements_array, string[] elements,
                                ref double[] a){
        /*
        Get the number of each element in a species from its lewis table data
        Inputs: 
             speciesList : array of strings of each species
             elements_array : array of associative arrays mapping species to elements
        Outputs:
            a : element_matrix mapping species to elemental composition 
        */
        ulong nsp,nel,j;
    
        nsp = speciesList.length;
        nel = elements.length; 
        a.length = nel*nsp;
        a[] = 0.0;
    
        foreach(i, atoms; enumerate(elements_array)){
            foreach (key, value; atoms){
                j = countUntil(elements, key);
                a[j*nsp + i] = value;
            }
        }
        return;
    }

    @nogc void pt_Assemble_Matrix(double[] a,double[] bi0,double[] G0_RTs,double p,double[] ns,double n,
                                  ref Matrix!double AB){
        /*
        Construct Iteration Matrix for reduced Newton Rhapson step, (eqn 2.24 and 2.26 from cea_I)
        */
        double lnns, lnn, lnp, akjaijnj, akjnjmuj, mus_RTj, bk;
        double nss, nsmus;
        int k,i,j,s;
        lnn = log(n);
        lnp = log(p/1e5); // Standard pressure for the tables is one BAR
    
        // Equation 2.24: k-> equation index, i-> variable index
        for (k=0; k<nel; k++){
            bk = 0.0; for (s=0; s<nsp; s++) bk += a[k*nsp + s]*ns[s];
            //A[k*neq + 0] = bk;
            AB[k,0] = bk;
    
            for (i=0; i<nel; i++){
                akjaijnj = 0.0;
                for (j=0; j<nsp; j++){
                    akjaijnj += a[k*nsp+j]*a[i*nsp+j]*ns[j];
                }
                //A[k*neq + i+1] = akjaijnj;
                AB[k,i+1] = akjaijnj;
            }
            akjnjmuj = 0.0;
            for (j=0; j<nsp; j++){
                if (ns[j]==0.0) continue;
                mus_RTj = G0_RTs[j] + log(ns[j]) - lnn + lnp;
                akjnjmuj += a[k*nsp+j]*ns[j]*mus_RTj;
    
            }
            //B[k] = bi0[k] - bk + akjnjmuj;
            AB[k,neq] = bi0[k] - bk + akjnjmuj;
    
            // Equation 2.26 - > (only the pii entries, we're highjacking k here to go across the last row)
            //A[nel*neq + k+1] = bk;
            AB[nel,k+1] = bk;
        }
        // Equation 2.26 - > (now the rest)
        nss = 0.0;
        nsmus = 0.0;
        for (j=0; j<nsp; j++){
            if (ns[j]==0.0) continue;
            mus_RTj = G0_RTs[j] + log(ns[j]) - lnn + lnp; // I guess its okay to compute this again
            nss += ns[j];
            nsmus += ns[j]*mus_RTj;
        }
        //A[nel*neq + 0]  = nss - n; // Yes, this is a problem. Shieeeeeet.
        //B[nel] = n - nss + nsmus;
        AB[nel,0]  = nss - n; 
        AB[nel,neq] = n - nss + nsmus;
    
        //for (i=0; i<neq; i++){
        //    for (j=0; j<neq; j++){
        //        printf("%f ", AB[i,j]);
        //    }
        //    printf("| %f\n", AB[i,neq]);
        //}
        return;
    }

    @nogc void pt_species_corrections(double[] S,double[] a,double[] G0_RTs,double p,double n,
                                       double[] ns, ref double[] dlnns){
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
            if (ns[s]==0.0) { dlnns[s] = 0.0; continue;}
            mu_RTs = G0_RTs[s] + log(ns[s]) - lnn + lnp;
    
            aispii = 0.0;
            for (i=0; i<nel; i++){
                aispii += a[i*nsp+s]*S[i+1]; // S[i+1] = pi_i, the lagrange multiplier
            }
            dlnns[s] = -mu_RTs + dlnn + aispii;
        }
        return; 
    }

    @nogc void pt_update_unknowns(double[] S,double[] dlnns,ref double[] ns, ref double n){
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
        double lnns,lnn,n_copy,lambda;
    
        lnn = log(n); 
        lambda = fmin(1.0, 0.5*fabs(lnn)/fabs(S[0]));
        n_copy = exp(lnn + lambda*S[0]); 
        n = n_copy;   // thing pointed to by n set to exp(lnn + S[0]);
    
        for (s=0; s<nsp; s++){
            if (ns[s]==0.0) continue;
            lnns = log(ns[s]);
    
            lambda = fmin(1.0, 0.5*fabs(lnn)/fabs(dlnns[s]));
            ns[s] = exp(lnns + lambda*dlnns[s]);
    
            if (ns[s]/n_copy<TRACELIMIT){
                ns[s] = 0.0;
                dlnns[s] = 0.0; // This species must be considered converged
            }
        }
    
        return;
    }
    
    @nogc int solve_pt(double p,double T,double[] X0,double[] M,double[] a,ref double[] X1,int verbose){
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
            verbose: print debugging information
    
        Output:
            X1 : Equilibrium Mole Fraction [nsp]  
        */
        int s,i,j,k,ntrace,errorcode;
        double M0,n,M1,errorL2,Rs;
    
        const double tol=1e-6;
        const int attempts=30;
    
        errorcode=0;
    
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
            Rs = Ru/M[s];
            auto c = _curves[s];
            G0_RTs[s] = c.eval_h(T)/Rs/T - c.eval_s(T)/Rs;
        }

        // Main Iteration Loop: 
        for (k=0; k<=attempts; k++){
            // 1: Perform an update of the equations
            pt_Assemble_Matrix(a, bi0, G0_RTs, p, ns, n, AB);
            try { gaussJordanElimination(AB); }
            catch (Exception caughtException) {
                for (s=0; s<nsp; s++) ns[s] = fmax(1e-3, ns[s]);
                continue;
            }
            
            for (j=0; j<neq; j++) S[j] = AB[j,neq];
            pt_species_corrections(S, a, G0_RTs, p, n, ns, dlnns);
            pt_update_unknowns(S, dlnns, ns, n);
    
            // Compute remaining error by checking species corrections
            errorL2 = 0.0;
            for (s=0; s<nsp; s++) errorL2 += dlnns[s]*dlnns[s];
            errorL2 = pow(errorL2, 0.5);
            if (verbose>0) printf("iter %d: %f %f %f %f  (%e) \n", k, n, ns[0], ns[1], ns[2], errorL2);
    
            // Check for imminent singularity
            ntrace=0;
            for (s=0; s<nsp; s++) if (ns[s]==0.0) ntrace++;
            if (ntrace==nsp-1) { // Pseudo convergence criteria, all the species but one are trace
                for (s=0; s<nsp; s++) if (ns[s]!=0.0) i=s;
                n = 1.0/M[i];
                ns[i] = n;
                errorL2 = 0.0;
                if (verbose>0) printf("Pseudo convergence! Remaining species: (%d)\n", i);
            }
    
            // Exit loop of convergence is achieved
            if (errorL2<tol) break;
    
            if (isNaN(errorL2)) {
                printf("Solver nan'd, exiting!\n");
                errorcode=1;
                break;
            }
    
            // Exit loop if too many attempts are undertaken
            if (k==attempts) {
                printf("Solver not converged, exiting!\n");
                errorcode=1;
                break;
            }
        }
        
        if ((verbose>0)&&(errorcode==0)) printf("Converged in %d iter, error: %e\n", k, errorL2);
        // Compute output composition
        M1 = 1.0/n;
        for (s=0; s<nsp; s++) X1[s] = M1*ns[s];
    
        // Don't need to clear memory manually
        return errorcode;
    }

    @nogc void rhou_Assemble_Matrix(double[] a,double[] bi0, double rho0, double u0, double T,
                                    double[] ns, double[] G0_RTs, double[] U0_RTs, double[] Cv0_Rs,
                                    ref Matrix!double AB){
        /*
        Construct Iteration Matrix for reduced Newton Rhapson rhou step, (eqn 2.45 and 2.47 from cea_I)
        */
        double mus_RTj, bk, u, akjaijnj, akjnjmuj, akjnjUj, njCvj, njUj2, njUjmuj, aijnjUj, Rs, Hs_RT;
        int k,i,j,s;
    
        u = 0.0;
        for (s=0; s<nsp; s++){
            Rs = Ru/M[s];
            auto c = _curves[s];

            // Note that Rowan uses Captial Cp,Cv to mean per mass, though I use it for per mole
            Hs_RT = c.eval_h(T)/Rs/T;
            G0_RTs[s] = Hs_RT - c.eval_s(T)/Rs;
            U0_RTs[s] = Hs_RT - 1.0;
            Cv0_Rs[s] = c.eval_Cp(T)/Rs - 1.0; 
            u += ns[s]*U0_RTs[s]*Ru*T;
        }
    
        // Equation 2.45: k-> equation index, i-> variable index
        for (k=0; k<nel; k++){
            bk = 0.0; for (s=0; s<nsp; s++) bk += a[k*nsp + s]*ns[s];
    
            for (i=0; i<nel; i++){
                akjaijnj = 0.0;
                for (j=0; j<nsp; j++){
                    akjaijnj += a[k*nsp+j]*a[i*nsp+j]*ns[j];
                }
                //A[k*neq + i+1] = akjaijnj; // Lagrange multiplier matrix entry
                AB[k,i+1] = akjaijnj; // Lagrange multiplier matrix entry
            }
            akjnjmuj = 0.0;
            akjnjUj = 0.0;
            for (j=0; j<nsp; j++){
                mus_RTj = G0_RTs[j] + log(rho0*ns[j]*Ru*T/1e5);
                akjnjmuj += a[k*nsp+j]*ns[j]*mus_RTj;
                akjnjUj  += a[k*nsp+j]*ns[j]*U0_RTs[j];
            }
            //A[k*neq + 0] = akjnjUj; // Temperature matrix entry
            //B[k] = bi0[k] - bk + akjnjmuj; // RHS of kth equation 2.45
            AB[k,0] = akjnjUj; // Temperature matrix entry
            AB[k,neq] = bi0[k] - bk + akjnjmuj; // RHS of kth equation 2.45
    
        }
        // Equation 2.47
        njCvj  = 0.0;
        njUj2  = 0.0;
        njUjmuj= 0.0;
    
        for (j=0; j<nsp; j++){
            mus_RTj = G0_RTs[j] + log(rho0*ns[j]*Ru*T/1e5);
            njCvj += ns[j]*Cv0_Rs[j];
            njUj2 += ns[j]*U0_RTs[j]*U0_RTs[j];
            njUjmuj += ns[j]*U0_RTs[j]*mus_RTj;
        }
        //A[nel*neq + 0]  = njCvj + njUj2;  // Temperature matrix entry
        //B[nel] = (u0 - u)/Ru/T + njUjmuj; // RHS 
        AB[nel,0]  = njCvj + njUj2;  // Temperature matrix entry
        AB[nel,neq] = (u0 - u)/Ru/T + njUjmuj; // RHS 
    
        for (i=0; i<nel; i++){
            aijnjUj = 0.0;
            for (j=0; j<nsp; j++){
                aijnjUj += a[i*nsp + j]*ns[j]*U0_RTs[j];
            }
            //A[nel*neq + 1+i] = aijnjUj; // Lagrange Multipliers       
            AB[nel,1+i] = aijnjUj; // Lagrange Multipliers       
        }
    
        //for (i=0; i<neq; i++){
        //    for (j=0; j<neq; j++){
        //        printf("%f ", AB[i,j]);
        //    }
        //    printf("| %f\n", AB[i,neq]);
        //}
        return;
    }

    @nogc void rhou_species_corrections(double[] S,double[] a,double[] G0_RTs,double[] U0_RTs,
                                        double rho0,double T, double[] ns, ref double[] dlnns){
        /*
        Compute delta_log(ns) from the reduced iteration equations from 
        equation 2.18m using the other deltas in S
        Inputs:
            S      : Corrections (pi1, pi2, pi3 ... dlog(n) [nel+1]
            a      : elemental composition array [nel,nsp]
            G0_RTs : Gibbs free energy of each species, divided by RT [nsp]
            U0_RTs : Internal energy of each species, divided by RT [nsp]
            rho0   : goal density (kg/m3)
            T      : current temperature guess (K)
            ns     : species moles/mixture kg [nsp]
    
        Outputs:
            dllns : change in log(ns) [nsp]
        */
        double dlnT,mus_RTs,aispii;
        int s,i;
        dlnT = S[0];
    
        // These should (should?) have been set during the assemble matrix step
        //for (s=0; s<nsp; s++){
        //    lp = lewis + 9*3*s;
        //    G0_RTs[s] = compute_G0_RT(T, lp);
        //    U0_RTs[s] = compute_H0_RT(T, lp) - 1.0;
        //}
        
        for (s=0; s<nsp; s++) {
            mus_RTs = G0_RTs[s] + log(rho0*ns[s]*Ru*T/1e5);
    
            aispii = 0.0;
            for (i=0; i<nel; i++){
                aispii += a[i*nsp+s]*S[i+1]; // S[i+1] = pi_i, the lagrange multiplier
            }
            dlnns[s] = -mus_RTs + U0_RTs[s]*dlnT + aispii;
        }
        return; 
    }

    @nogc void rhou_update_unknowns(double[] S,double[] dlnns, ref double[] ns, ref double T){
        /*
        Add corrections to unknown values (ignoring lagrange multipliers)
        Inputs:
            S : vector of corrections from matrix step [nel+1]
            dlnns : vector of species mole/mixture corrections [nsp]
        Outputs:
            ns : vector of species mole/mixtures [nsp]
            T  : pointer to temperature guess (passed by reference!) [1]
        */
        int s;
        double lnns,lnT,n,lnn,lambda;
        lnT = log(T);  // Removed dereference operation
        lambda = fmin(1.0, 0.5*fabs(lnT)/fabs(S[0]));
        T = exp(lnT + lambda*S[0]); // removed dereference operation 
    
        n = 0.0; for (s=0; s<nsp; s++) n+=ns[s]; lnn=log(n);
    
        for (s=0; s<nsp; s++){
            lnns = log(ns[s]);
            lambda = fmin(1.0, fabs(lnn)/fabs(dlnns[s]));
            ns[s] = exp(lnns + lambda*dlnns[s]);
    
            if (ns[s]/n<TRACELIMIT){
                ns[s] = 0.0;
                dlnns[s] = 0.0; // This species must be considered converged
            }
        }
        return;
    }

    @nogc int solve_rhou(double rho,double u,double[] X0,double[] M,double[] a,
                         ref double[] X1, ref double Teq, int verbose){
        /*
        Compute the equilibrium composition X1 at a fixed volume and internal energy 
        Inputs:
            rho   : target Density (kg/m3)
            u     : target internal energy (J/kg)
            X0    : Intial Mole fractions [nsp]
            M     : Molar Mass of each species (kg/mol) [nsp]
            a     : elemental composition array [nel,nsp]
            verbose: print debugging information
    
        Output:
            X1  : Equilibrium Mole Fraction [nsp]  
            Teq: Equilibrium Temperature 
        */
        int s,i,j,k,errorcode;
        double M0,n,M1,errorL2,T;
    
        const double tol=1e-6; //FIXME: Make this a class variable?
        const int attempts=10; //FIXME: Make this a class variable?
    
        errorcode=0;
    
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
    
        //T = temperature_guess(nsp, u, M0, X0, lewis); // Make sure this is set on startup
        T = Teq; // Should copy? No?
        if (verbose>0) printf("Guess T: %f\n", T);
    
        // Begin Iterations
        for (k=0; k<attempts; k++){
            rhou_Assemble_Matrix(a,bi0,rho,u,T,ns,G0_RTs,U0_RTs,Cv0_Rs,AB);
            try { gaussJordanElimination(AB); }
            catch (Exception caughtException) {
                for (s=0; s<nsp; s++) ns[s] = fmax(1e-3, ns[s]);
                continue;
            }

            for (j=0; j<neq; j++) S[j] = AB[j,neq];
            rhou_species_corrections(S,a,G0_RTs,U0_RTs,rho,T,ns,dlnns);
            rhou_update_unknowns(S, dlnns, ns, T);
    
            errorL2 = 0.0;
            for (s=0; s<nsp; s++) errorL2 += dlnns[s]*dlnns[s];
            errorL2 = pow(errorL2, 0.5);
            if (verbose>0) printf("iter %d: %f %f %f %f  (%e) \n", k, T, ns[0], ns[1], ns[2], errorL2);
            if (errorL2<tol) break;
    
            if (k>=attempts) {
                printf("Solver not converged, exiting!\n");
                errorcode=1;
                break;
            }
        }
        
        if ((verbose>0)&&(errorcode==0)) printf("Converged in %d iter, error: %e\n", k, errorL2);
        // Compute output composition
        n = 0.0;
        for (s=0; s<nsp; s++) n += ns[s];
        M1 = 1.0/n;
        for (s=0; s<nsp; s++) X1[s] = M1*ns[s];
        Teq = T;
    
        // Don't need to clear memory manually
        return errorcode;
    }

} // end class EQGas 


version(eq_gas_test) {
    int main() {
        string[3] species = ["CO2", "CO", "O2"];
        double[3] Y0 = [1.0, 0.0, 0.0];
        double[3] Xt = [0.66010397, 0.22659735, 0.11329868];
        double[3] Xf = [0.0, 0.0, 0.0];
        auto gm = new EqGas("sample-data/co2_test.lua");
        auto Q = new GasState(3, 0);
        auto Q2= new GasState(3, 0);

        Q.massf = Y0;
        Q.p = 0.1*101.35e3;
        Q.T = 2500.0;
        gm.update_thermo_from_pT(Q);
        gm.massf2molef(Q, Xf); 
        writeln("Test pt: ");
        writeln("    Found: ", Xf);
        writeln("    Target: ", Xt);
        writeln("    ", Q);
        writeln("    ");

        Q2.massf = Y0;
        Q2.rho = Q.rho;
        Q2.u = Q.u;
        Q2.T = Q.T*1.2;
        gm.update_thermo_from_rhou(Q2);
        gm.massf2molef(Q2, Xf); 
        writeln("Test rhou: ");
        writeln("    Found: ", Xf);
        writeln("    Target: ", Xt);
        writeln("    ", Q2);

        // Test perfectly converged case
        auto Q3= new GasState(3, 0);
        Q3.massf[0] = 1.0;
        Q3.massf[1] = 0.0;
        Q3.massf[2] = 0.0;
        Q3.rho = 0.107294;
        Q3.u = -8999109.972193;
        Q3.T = 296.0;
        Q3.p = 6000.0;
        gm.update_thermo_from_rhou(Q3);
        writeln("Test converged rhou: ");
        writeln(Q3);

        return 0;
    } // end main()
}
