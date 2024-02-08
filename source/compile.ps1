gcc -I. -fPIC -Wall -std=c99 -O3 -c .\thermo.c
gcc -I. -fPIC -Wall -std=c99 -O3 -c .\linalg.c
gcc -I. -fPIC -Wall -std=c99 -O3 -c .\common.c
gcc -I. -fPIC -Wall -std=c99 -O3 -c .\pt.c
gcc -I. -fPIC -Wall -std=c99 -O3 -c .\rhot.c
gcc -I. -fPIC -Wall -std=c99 -O3 -c .\rhou.c
gcc -I. -fPIC -Wall -std=c99 -O3 -c .\ps.c
gcc -I. -fPIC -Wall -std=c99 -O3 -c .\eqc.c

gcc -I. -fPIC -Wall -std=c99 -O3 -shared .\thermo.o .\linalg.o .\common.o .\pt.o .\rhou.c .\rhot.o .\ps.o .\eqc.o -lm -o .\libeqc.dll
