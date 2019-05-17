/*
Test module for d integration into eilmer 4

@author: Nick Gibbons
*/

import std.stdio;
import std.string;
import std.file;
import std.conv;
import std.algorithm: canFind, countUntil;
import std.algorithm.sorting: sort;
import std.range: enumerate, take;
import ceq;

const string DBPATH = "./test.inp";

int readdb(string species, ref string[] lewisdata){
    /* 
    Read a species raw data from the lewis thermodynamic database.

    Inputs: 
        species : species formula
    Outputs:
        lewisdata : array of strings representing each line of the database entry (should start blank!)
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
        M : molecular mass in kg/mol (that's gmols, unlike in the lewis database, which uses kmols)
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

void startup(string[] speciesList,
             ref string[] element_set, ref double[] lewiscoeffs, ref double[string][] element_map,
             ref double[] a, ref double[] M){
    /*
    Entry point function for equilibrium chemistry simulations
    Inputs: 
         speciesList : array of strings of each species
    Outputs:
        element_set : array of elements in the complete chemical system
        lewiscoeffs : numbers for computing thermodynamic data from lewis tables
        element_map : array of associative arrays mapping species to elements
        a : element_matrix mapping species to elemental composition 
        M : array of molecular weights
    */
    string[] species_lewisdata;
    double[] species_lewiscoeffs;
    double[string] species_elements;

    foreach(species; speciesList){
        species_lewisdata.length=0;
        species_lewiscoeffs.length=0;
        species_elements.clear();

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
    return;
}

int main(string[] args){
    int nsp,nel,code;
    double[] Xs0;
    double[] Xst;
    double[] Xs1;
    double[] M;
    string[] speciesList;
    string[] element_set;
    double[] lewiscoeffs;
    double[] a;
    double[string][] element_map;

    double u1,T1,p1;

    nsp = 3;
    speciesList.length = nsp;
    speciesList[0] = "CO2";
    speciesList[1] = "CO";
    speciesList[2] = "O2";

    Xs0.length = nsp; 
    Xs0[0] = 1.0;
    Xs0[1] = 0.0;
    Xs0[2] = 0.0;

    Xst.length = nsp; 
    Xst[0] = 0.66108962603325838;
    Xst[1] = 0.22594024931116111;
    Xst[2] = 0.11297012465558055;

    Xs1.length = nsp;
    Xs1[] = 0.0;

    startup(speciesList, element_set, lewiscoeffs, element_map, a, M);
    nel = to!int(element_set.length);

    writeln("elements", element_set);
    writeln("lewiscoeffs", lewiscoeffs);
    writeln("element_map", element_map);
    writeln("a", a);
    writeln("M", M);

    T1 = 2500.0;
    u1 = 0.0;
    p1= 0.1*101.35e3;

    u1 = get_u(T1, Xst.ptr, nsp, lewiscoeffs.ptr, M.ptr);
    writeln("testu: ", u1);

    code = pt(p1,T1,Xs0.ptr,nsp,nel,lewiscoeffs.ptr,M.ptr,a.ptr,Xs1.ptr,1);
    writeln("Done pt: ", code);
    writeln("Xs1: ", Xs1);
    writeln("Xst: ", Xst);

    return 0;
}
