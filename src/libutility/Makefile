#**
# Makefile for libutility.a, a library for all utility things 
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

OBJS = 	alloc_struct.o \
       	cosmology.o\
	general.o\
	loadbalance.o\
	specific.o \
    util_bsearch.o \
    timer.o \
	darkenergy.o

# Catch empty make call
all:
	${MAKE} libutility.a

# Cleaning up
clean:
	rm -f libutility.a $(OBJS)

# The library itself
libutility.a: $(OBJS)
	${AR} -r libutility.a $(OBJS)


# The individual routines stuff
alloc_struct.o: alloc_struct.c
	$(CC) $(CFLAGS) -c -o alloc_struct.o alloc_struct.c

cosmology.o: cosmology.c
	$(CC) $(CFLAGS) -c -o cosmology.o cosmology.c

darkenergy.o: darkenergy.c
	$(CC) $(CFLAGS) -c -o darkenergy.o darkenergy.c

general.o: general.c
	$(CC) $(CFLAGS) -c -o general.o general.c

loadbalance.o: loadbalance.c
	$(CC) $(CFLAGS) -c -o loadbalance.o loadbalance.c

specific.o: specific.c
	$(CC) $(CFLAGS) -c -o specific.o specific.c

util_bsearch.o: util_bsearch.c util_bsearch.h
	$(CC) $(CFLAGS) -c -o util_bsearch.o util_bsearch.c

timer.o: timer.c timer.h
	$(CC) $(CFLAGS) -c -o timer.o timer.c
