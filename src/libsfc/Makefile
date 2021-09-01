#**
# Makefile for libsfc.a, a small library providing space filling curve
# algorithms.
#
# version 0.2 (01.12.2006), Steffen Knollmann (sknollmann@aip.de)
#
# Hilbert Curve implementation copyright 1998, Rice University
#**

# If ${OPTIMIZE}, ${CCFLAGS} and/or ${DEFINEFLAGS} are set, they will be
# used for the compilation
CFLAGS=${OPTIMIZE} ${CCFLAGS} ${DEFINEFLAGS}

# If a different MAKE ist specified use that one, otherwise try you luck
# with just make; maybe it will work.
MAKE?=make

# The same for ar
AR?=ar

# PHONY targets
.PHONY: all clean

# Catch empty make call
all:
	${MAKE} libsfc.a

# Cleaning up
clean:
	rm -f libsfc.a sfc_curve.o sfc_boundary.o hilbert_util.o hilbert.o

# The library itself
libsfc.a: sfc_curve.o sfc_boundary.o hilbert_util.o hilbert.o
	${AR} -r libsfc.a sfc_curve.o sfc_boundary.o hilbert_util.o hilbert.o

# A generic interface to space filling curves
sfc_curve.o: hilbert_util.h sfc_curve.h sfc_curve.c
	$(CC) $(CFLAGS) -c -o sfc_curve.o sfc_curve.c

# Used to construct the boundary of space fillung curves
sfc_boundary.o: sfc_curve.h sfc_boundary.h sfc_boundary.c
	$(CC) $(CFLAGS) -c -o sfc_boundary.o sfc_boundary.c

# A utility modul for Hilbert Curves
hilbert_util.o: hilbert.h hilbert_util.h hilbert_util.c
	$(CC) $(CFLAGS) -c -o hilbert_util.o hilbert_util.c

# The Hilbert Curve implementation by Doug Moore
hilbert.o: hilbert.h hilbert.c
	$(CC) $(CFLAGS) -c -o hilbert.o hilbert.c
