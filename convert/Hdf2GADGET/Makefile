
include make.COSMA4

OPTIM = -O2 -Mvect
F90FLAGS = $(DEBUG) #-convert big_endian -O #-fast
FFLAGS = $(FDEBUG)
CFLAGS = $(CDEBUG)


.SUFFIXES:
.SUFFIXES:  .exe .o .for .f .f90

.for.o:
	$(FC) $(INC) -c $< -o $*.o $(F77FLAGS)  

.f.o:
	$(FC) $(INC) -c $< -o $*.o $(F77FLAGS)  

.o.exe:
	$(FC) $(INC) $< -o $*.exe $(F90FLAGS)   

.f90.o:
	$(FC) $(INC) -c $< -o $*.o $(F90FLAGS)

Hdf2GADGET: Hdf2GADGET_io.o Hdf2GADGET_subroutines.o Hdf2GADGET.o
	$(FC) $(INC) $(FFLAGS) $(OPTIONS) $(LIB) Hdf2GADGET_io.o Hdf2GADGET_subroutines.o Hdf2GADGET.o -o Hdf2GADGET $(F90FLAGS)

clean:
	rm -rf *.o *.mod *exe



