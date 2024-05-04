CC=g++
CFLAGS=-Iutilities/ -Iinclude -fopenmp
LFLAGS=-larmadillo -fopenmp
ODIR=lib

all: main

main: run/main.cpp $(ODIR)/CNDO2.o $(ODIR)/timer.o
	$(CC) -o $@ $^ $(CFLAGS) $(LFLAGS)

$(ODIR)/CNDO2.o: source/CNDO2.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/timer.o: utilities/timer.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o main