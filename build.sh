#! /bin/bash
gcc -c -fpic thermo.c 
gcc -c -fpic linalg.c 
gcc -c -fpic pt.c 
gcc -c -fpic rhou.c 
gcc -c -fpic ceq.c
gcc -shared thermo.o linalg.o pt.o rhou.o ceq.o -lm -o libceq.so
dmd deq.d thermo.o linalg.o pt.o rhou.o ceq.o
rm *.o 
