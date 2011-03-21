CC = g++
#CFLAGS = -O3

TARGETS = orig serial-hfk

all:	$(TARGETS)

serial-hfk: serial-hfk.cpp
	$(CC) -o $@ $(LIBS) serial-hfk.cpp

orig: orig.cpp
	$(CC) -o $@ orig.cpp

clean:
	rm -f *.o $(TARGETS)
