#! /bin/bash
gcc -c -fpic thermo.c 
gcc -c -fpic linalg.c 
gcc -c -fpic ceq.c
gcc -shared thermo.o linalg.o ceq.o -lm -o libceq.so
rm *.o 
