CC = g++
LINBOX =  -L/usr/local/lib -L/usr/local/lib64 -llinbox -lgivaro -lgmpxx -lgmp -lcblas -llapack
#CFLAGS = -O3

TARGETS = orig serial-hfk low-mem sparse-linbox serial-hfk-split

all:	$(TARGETS)

sparse-linbox: sparse-linbox.cpp
	$(CC) -o $@ $(LIBS) $(LINBOX) sparse-linbox.cpp

low-mem: low-mem.cpp matrix-rep.h
	$(CC) -o $@ $(LIBS) $(LINBOX) low-mem.cpp

simple-test: simple-test.cpp simple-mat.h
	$(CC) -o $@ $(LIBS) $(LINBOX) simple-test.cpp

serial-hfk-split: serial-hfk-split.cpp
	$(CC) -o $@ $(LIBS) serial-hfk-split.cpp

serial-hfk: serial-hfk.cpp
	$(CC) -o $@ $(LIBS) serial-hfk.cpp

orig: orig.cpp
	$(CC) -o $@ orig.cpp

clean:
	rm -f *.o $(TARGETS)
