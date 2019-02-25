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
import std.range: enumerate;

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
    i=0;
    foreach(line; fiter){
       if (i==nlines) break;
       lewisdata ~= to!string(nextline);
       i++;
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
    M /= 1000.0;
    return M;
}

void species_elements(string[] lewisdata, ref double[string] elements){
    /*
    Get the number of each element in a species from its lewis table data
    Inputs: 
        lewisdata : array of strings representing each line of the database entry 
    Outputs:
        elements : associative array (dict) of element types (should be blank)
    */
    if (elements.length!=0) throw new Exception("elements.length!=0");
    int i,n,s,e;
    double nel;
    string element,entry;
    n = 1;
    
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
    Get the number of each element in a species from its lewis table data
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

// TODO: Make a separate function to get "elements" and sort it
void compute_element_matrix(string[] speciesList, double[string][] elements_array, ref double[] a){
    /*
    Get the number of each element in a species from its lewis table data
    Inputs: 
         speciesList : array of strings of each species
         elements_array : array of associative arrays mapping species to elements
    Outputs:
        a : element_matrix mapping species to elemental composition (Return the element list too?)
    */
    ulong nsp,nel,j;
    string[] elements;
    string el;

    foreach(e; elements_array){
        foreach(key; e.keys()){
            if (!elements.canFind(key)) elements ~= key;
        }
    }
    elements.sort();

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

int main(string[] args){
    int nsp,i,s,e;
    double[] Xs0;
    double[] Xst;
    double[] M;
    string[] speciesList;
    string[] lewisdata;
    double[] lewiscoeffs;
    double[] lewiscoeffs_array;
    double[] a;
    double[string] elements;
    double[string][] elements_array;

    nsp = 3;
    speciesList.length = nsp;
    speciesList[0] = "CO2";
    speciesList[1] = "CO";
    speciesList[2] = "O2";

    Xs0.length = nsp; 
    Xs0[0] = 1.0;
    Xs0[1] = 0.0;
    Xs0[0] = 0.0;

    Xst.length = nsp; 
    Xst[0] = 0.66108962603325838;
    Xst[1] = 0.22594024931116111;
    Xst[0] = 0.11297012465558055;

    for (i=0; i<nsp; i++){
        lewisdata.length=0;
        lewiscoeffs.length=0;
        elements.clear();

        readdb(speciesList[i], lewisdata);

        M ~= species_molecular_mass(lewisdata); 
        species_elements(lewisdata, elements);
        elements_array ~= elements.dup();
        species_thermo_coefficients(lewisdata, lewiscoeffs);
        if (lewiscoeffs.length==18) lewiscoeffs ~= lewiscoeffs[9..18]; //safety catch for nlines=2
        lewiscoeffs_array ~= lewiscoeffs.dup();
    }
    compute_element_matrix(speciesList, elements_array, a);
    writeln("a", a);

    return 0;
}
