CC=gcc
DMD=dmd
CFLAGS=-I. -fpic

libceq: thermo.o linalg.o pt.o rhou.o ceq.o
	$(CC) -shared thermo.o linalg.o pt.o rhou.o ceq.o -lm -o libceq.so

deq: thermo.o linalg.o pt.o rhou.o ceq.o
	$(DMD) deq.d thermo.o linalg.o pt.o rhou.o ceq.o

clean: 
	rm *.o 
