FC = gfortran
CC = gcc
FFLAGS = 
CFLAGS = -g -fpic
LINKFLAGS = -I/Users/aknebe/bin/hdf5//lib -I/Users/aknebe/bin/hdf5//include/ -I/Users/aknebe/Office/Source/HaloFinder/ahf/convert/Hdf2GADGET/hdf5_wrapper/src/../lib -L/Users/aknebe/Office/Source/HaloFinder/ahf/convert/Hdf2GADGET/hdf5_wrapper/src/../lib -I/Users/aknebe/Office/Source/HaloFinder/ahf/convert/Hdf2GADGET/hdf5_wrapper/src/../lib -lhdfwrapper

.SUFFIXES:
.SUFFIXES: .f .c .o .f90

ALLFFLAGS = $(FFLAGS) $(LINKFLAGS)
ALLCFLAGS = $(CFLAGS)

.f.o:
	$(FC) $(ALLFFLAGS) -c $< -o $*.o
.f90.o:
	$(FC) $(ALLFFLAGS) -c $< -o $*.o
.c.o:
	$(CC) $(ALLCFLAGS) -c $< -o $*.o

all: test_write.exe test_read.exe test_extend.exe

test_write.exe: test_write.o
	$(FC)   test_write.o $(ALLFFLAGS) -o test_write.exe

test_read.exe: test_read.o
	$(FC)  test_read.o $(ALLFFLAGS) -o test_read.exe

test_extend.exe: test_extend.o
	$(FC)  test_extend.o $(ALLFFLAGS) -o test_extend.exe

clean:
	rm -f *.o *.mod
	rm -f test_read.exe test_write.exe test_extend.exe


