#! /bin/bash
gcc -c thermo.c 
gcc -c linalg.c 
gcc -c ceq.c
gcc thermo.o linalg.o ceq.o -lm -o ceq
rm *.o 
