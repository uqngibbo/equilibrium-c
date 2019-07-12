CC ?= gcc
CFLAGS := -I. -fPIC
INSTALL_DIR=/home/qungibbo/programs/ceq

thermo.o: thermo.c thermo.h
	$(CC) $(CFLAGS) -c thermo.c

linalg.o: linalg.c linalg.h
	$(CC) $(CFLAGS) -c linalg.c

common.o: common.c common.h
	$(CC) $(CFLAGS) -c common.c

rhou.o: rhou.c thermo.h linalg.h common.h rhou.h
	$(CC) $(CFLAGS) -c rhou.c

pt.o: pt.c thermo.h linalg.h common.h pt.h
	$(CC) $(CFLAGS) -c pt.c

ceq.o: ceq.c thermo.h rhou.h pt.h ceq.h
	$(CC) $(CFLAGS) -c ceq.c


ceq: thermo.o linalg.o common.o pt.o rhou.o ceq.o
	$(CC) $(CFLAGS) -shared thermo.o linalg.o common.o pt.o rhou.o ceq.o -lm -o libceq.so

install: 
	cp libceq.so $(INSTALL_DIR)
	sed "s+DBPATH='../thermo.inp'+DBPATH='$(INSTALL_DIR)/thermo.inp'+" pyeq.py > pyeq1.py
	sed "s+LIBPATH='./libceq.so'+LIBPATH='$(INSTALL_DIR)/libceq.so'+" pyeq1.py > pyeq2.py
	cp pyeq2.py $(INSTALL_DIR)/pyeq.py
	rm -rf pyeq1.py pyeq2.py
	cp '../thermo.inp' $(INSTALL_DIR)
	rm -rf $(INSTALL_DIR)/tests
	cp -r ../tests $(INSTALL_DIR)

moduletests: thermo.o linalg.o
	$(CC) thermo.c -D TEST -lm -o ../tests/testthermo
	$(CC) linalg.c -D TEST -lm -o ../tests/testlinalg

clean: 
	rm -rf *.o 
	rm -rf libceq.so
