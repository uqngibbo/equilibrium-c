#! /bin/bash
gcc thermo.c -D TEST -lm -o testthermo
gcc linalg.c -D TEST -lm -o testlinalg
