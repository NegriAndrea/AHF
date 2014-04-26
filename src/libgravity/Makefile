#**
# Makefile for libgravity.a, a library for all things related to the gravity solver
#
# version 0.0 (09.09.2008), Alexander Knebe (aknebe@aip.de)
#**

# If ${OPTIMIZE}, ${CCFLAGS} and/or ${DEFINEFLAGS} are set, they will be
# used for the compilation
CFLAGS+=${OPTIMIZE} ${CCFLAGS} ${DEFINEFLAGS}

# If a different MAKE is specified use that one, otherwise try your luck
# with just make; maybe it will work.
MAKE?=make

# The same for ar
AR?=ar

# PHONY targets
.PHONY: all clean

OBJS = 	fft_potential.o \
			gs.o\
			residual.o\
			solve_gravity.o\
			trunc_err.o

# Catch empty make call
all:
	${MAKE} libgravity.a

# Cleaning up
clean:
	rm -f libgravity.a $(OBJS)

# The library itself
libgravity.a: $(OBJS)
	${AR} -r libgravity.a $(OBJS)


# The individual routines stuff
fft_potential.o: fft_potential.c
	$(CC) $(CFLAGS) -c -o fft_potential.o fft_potential.c

gs.o: gs.c
	$(CC) $(CFLAGS) -c -o gs.o gs.c

residual.o: residual.c
	$(CC) $(CFLAGS) -c -o residual.o residual.c

solve_gravity.o: solve_gravity.c
	$(CC) $(CFLAGS) -c -o solve_gravity.o solve_gravity.c

trunc_err.o: trunc_err.c
	$(CC) $(CFLAGS) -c -o trunc_err.o trunc_err.c


