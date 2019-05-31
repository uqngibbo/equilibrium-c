/**
 * eq_gas.d
 * Equilibrium gas mixture based on ceq library
 *
 * Author: Nick Gibbons 
 * Notes: - Start
 */

module gas.eq_gas;

import std.math;
import std.stdio;
import std.string;
import std.conv : to;
import util.lua;
import util.lua_service;
import std.file;
import std.algorithm: canFind, countUntil;
import std.algorithm.sorting: sort;
import std.range: enumerate, take;
import nm.number;
import nm.bbla;
import gas.thermo.cea_thermo_curves;
import gas.gas_model;
import gas.gas_state;

const string DBPATH = "/home/qungibbo/programs/ceq/thermo.inp";
const double Ru=8.3144621;      // Universal Gas Constant (J/mol/K)
const double TRACELIMIT=1e-7;   // Trace species limiter (for ns/n)

class EqGas: GasModel {
public:
    this(string mixtureName, string[] speciesList){
        version(complex_numbers) {
            throw new Error("Do not use with complex numbers.");
        }

        // We're going to trick the flow solver into thinking that we have only one species
        _n_modes = 0;
        _n_species = 1;
        _species_names.length =  _n_species;
        _species_names[0] = mixtureName;

        // Now set up our own secret variables
        string[] species_lewisdata;
        double[] species_lewiscoeffs;
        double[string] species_elements;

        foreach(species; speciesList){
            species_lewisdata.length=0;
            species_lewiscoeffs.length=0;
            species_elements.clear();

            species_list ~= species;
            readdb(species, species_lewisdata);

            M ~= species_molecular_mass(species_lewisdata); 

            species_elements_from_lewisdata(species_lewisdata, species_elements);
            element_map ~= species_elements.dup();

            species_thermo_coefficients(species_lewisdata, species_lewiscoeffs);
            if (species_lewiscoeffs.length==18) species_lewiscoeffs ~= species_lewiscoeffs[9..18]; 
            lewiscoeffs ~= species_lewiscoeffs.dup();
        }

        compile_element_set(element_map, element_set);
        compile_element_matrix(speciesList, element_map, element_set, a);

        nel = to!int(element_set.length);
        nsp = to!int(speciesList.length);
        _mol_masses.length = nsp;
        _mol_masses = M;  // FIXME: this could be a problem, put the calc in the getter ??
        X0.length = nsp; X1.length = nsp;

        neq= nel+1;
        S.length      = neq;     // Iteration unknown vector
        G0_RTs.length = nsp;     // Species Gibbs Free Energy
        ns.length     = nsp;     // Species moles/mixture mass
        bi0.length    = nel;     // starting composition coefficients
        dlnns.length  = nsp;     // starting composition coefficients
        AB = new Matrix!double(neq, neq+1);
        
        // FIXME: eventually this will be moved to the other constructor and not have a hardcoded name
        auto L = init_lua_State();
        doLuaFile(L, "sample-data/co2_test.lua");
        foreach ( isp; 0..nsp ) {
            lua_getglobal(L, "db");
            lua_getfield(L, -1, speciesList[isp].toStringz);
            lua_getfield(L, -1, "thermoCoeffs");
            _curves ~= createCEAThermo(L, Ru/M[isp]);
            lua_pop(L, 1);
            lua_pop(L, 1);
            lua_pop(L, 1);
        }

        // TODO you need viscous properties setup here. 
    } // End Basic contructor


    this(lua_State* L) {
        lua_getglobal(L, "EqGas");

        // getXXX functions automatically clear the LUA stack, reseting it to EQGas
        string mixtureName = getString(L, -1, "mixtureName");

        string[] speciesList;
        getArrayOfStrings(L, -1, "speciesList", speciesList);

        this(mixtureName, speciesList);
    } // end constructor from a Lua file


    override string toString() const
    {
        return "EqGas(species=[TODO])";
    }

    override void update_thermo_from_pT(GasState Q) 
    {
        int error;

        massf2molef(Q, X0); 
        error = solve_pt(Q.p,Q.T,X0,M,a,X1,1);
        molef2massf(X1, Q);
        // Set everything else in Q that isn't p,t 

    }
    override void update_thermo_from_rhou(GasState Q)
    {
        //int error;
        //double* Tp;
        //double Mmix;
        //massf2molef(Q, X0); 
        //error = solve_rhou(Q.rho,Q.u,X0.ptr,nsp,nel,lewiscoeffs.ptr,M.ptr,a.ptr,X1.ptr,Tp,0);

        //// Set everything else in Q that isn't rho,u 
        //molef2massf(X1, Q);
        //Mmix = 0.0; for(int s; s<nsp; s++) Mmix += X1[s]*M[s];

        //Q.T = *Tp;
        //Q.p = Ru/Mmix*Q.rho*Q.T;
        //Q.p_e = Q.p;
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
    override void update_sound_speed(GasState Q) const
    {
        // FIXME

    }
    override void update_trans_coeffs(GasState Q)
    {
        // FIXME
    }

    override number dudT_const_v(in GasState Q) const
    {
        return to!number(Q.ceaSavedData.Cp - Q.ceaSavedData.Rgas);
    }
    override number dhdT_const_p(in GasState Q) const
    {
        return to!number(Q.ceaSavedData.Cp);
    }
    override number dpdrho_const_T(in GasState Q) const
    {
        return Q.ceaSavedData.Rgas * Q.T;
    }
    override number gas_constant(in GasState Q) const
    {
        return to!number(Q.ceaSavedData.Rgas);
    }
    override number internal_energy(in GasState Q) const
    {
        return Q.u;
    }
    override number enthalpy(in GasState Q) const
    {
        return to!number(Q.ceaSavedData.h);
    }
    override number entropy(in GasState Q) const
    {
        return to!number(Q.ceaSavedData.s);
    }

    @nogc @property override uint gasstate_n_species() const { return nsp; }

private:
    int nel, nsp, neq;
    string[] species_list;
    string[] element_set;
    double[] lewiscoeffs;
    double[string][] element_map;
    double[] a;
    double[] M;
    double[] X0,X1;
    double[] S, G0_RTs, ns, bi0, dlnns; // Dynamic arrays
    Matrix!double AB;
    CEAThermo[] _curves;

    int readdb(string species, ref string[] lewisdata){
        /* 
        Read a species raw data from the lewis thermodynamic database.
    
        Inputs: 
            species : species formula
        Outputs:
            lewisdata : array of strings representing each line of the database entry 
        */ 
        int i,nlines,nentries;
        string sp_with_space = species ~ " ";
        auto f = File(DBPATH);
        auto fiter = f.byLine();
        char[] nextline;
        
        if (lewisdata.length!=0) throw new Exception("lewisdata.length!=0");
    
        // Search through the file looking for 'species'
        foreach(line; fiter) if (line.startsWith(sp_with_space)) break;
    
        // Check to see if we found it
        if (f.eof()) throw new Exception("Species " ~ species ~ " not found!");
    
        // If we've arrived here we're just before the database entry, check how many entries
        fiter.popFront();
        nextline = fiter.front;
        nentries = to!int(nextline.split()[0]);
        nlines = nentries*3+1;
    
        // Now compilte them into a dynamic array, including the "nextline" from above
        foreach(line; fiter.take(nlines)){
           lewisdata ~= to!string(line);
        }
        return 0;
    }
    
    double species_molecular_mass(string[] lewisdata){
        /*
        Get the molecular weight of a species from its lewis table data.
        Inputs: 
            lewisdata : array of strings representing each line of the database entry 
        Outputs:
            M : molecular mass in kg/mol (that's kmols, unlike in the lewis database, which uses gmols)
        */
        double M;
    
        auto words = lewisdata[0].split();
        M = to!double(words[$-2]);
        M /= 1000.0; // Change to gram-moles, (or regular moles...)
        return M;
    }
    
    void species_elements_from_lewisdata(string[] lewisdata, ref double[string] elements){
        /*
        Get the number of each element in a species from its lewis table data
        Inputs: 
            lewisdata : array of strings representing each line of the database entry 
        Outputs:
            elements : associative array (dict) of element types (should start empty)
        */
        if (elements.length!=0) throw new Exception("elements.length!=0");
        int i,s,e;
        double nel;
        string element,entry;
        
        for (i=0; i<5; i++){
            s = 10+8*i;
            e = 10+8*(i+1);
            entry = lewisdata[0][s..e];
    
            auto entrysplit = entry.split();
            if (entrysplit.length!=2) continue;
    
            element = to!string(entrysplit[0]);
            nel = to!double(entrysplit[1]);
    
            elements[element] = nel;
        }
        return;
    }
    
    void species_thermo_coefficients(string[] lewisdata, ref double[] thermo_coefficients){
        /*
        Get the numbers representing thermodynamic property curve fits from the lewis table data
        Inputs: 
            lewisdata : array of strings representing each line of the database entry 
        Outputs:
            thermo_coefficients : numbers for computing thermodynamic data from lewis tables
        */
        if (thermo_coefficients.length!=0) throw new Exception("thermo_coefficients.length!=0");
        int i,j,s,e,nentries,nlines;
        string line,D,E;
    
        nlines = to!int(lewisdata.length);
        nentries = (nlines-1)/3;
    
        for (i=0; i<nentries; i++){
            line = lewisdata[2+3*i].replace("D","E"); // Fix Fortran Double nonsense
            for (j=0; j<5; j++){
                s = 16*j;
                e = 16*(j+1);
                thermo_coefficients ~= to!double(line[s..e].strip());
            }
    
            line = lewisdata[2+3*i+1].replace("D","E"); // Fix Fortran Double nonsense
            for (j=0; j<5; j++){
                if (j==2) continue; // Skip blank entry
                s = 16*j;
                e = 16*(j+1);
                thermo_coefficients ~= to!double(line[s..e].strip());
            }
        }
        return;
    }
    
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
    
            lambda = fmin(1.0, fabs(lnn)/fabs(dlnns[s]));
            ns[s] = exp(lnns + lambda*dlnns[s]);
    
            if (ns[s]/n_copy<TRACELIMIT) ns[s] = 0.0;
        }
    
        return;
    }
    
    // TODO: Change io variables from mole to mass fractions maybe?
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
        double M0,n,M1,errorL2,thing,Rs;
    
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
            gaussJordanElimination(AB);
            
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

} // end class EQGas 


version(eq_gas_test) {
    int main() {
        string[3] species = ["CO2", "CO", "O2"];
        double[3] Y0 = [1.0, 0.0, 0.0];
        double[3] Xt = [0.66010397, 0.22659735, 0.11329868];
        double[3] Xf = [0.0, 0.0, 0.0];
        auto gm = new EqGas("co2_eq", species);
        auto Q = new GasState(3, 0);
        Q.massf = Y0;
        Q.p = 0.1*101.35e3;
        Q.T = 2500.0;
        writeln(Q);
        gm.update_thermo_from_pT(Q);
        gm.massf2molef(Q, Xf); 
        writeln("Found: ", Xf);
        writeln("Target: ", Xt);

        return 0;
    } // end main()
}
