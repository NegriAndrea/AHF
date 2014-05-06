#**
# Makefile for libahf.a, a library for all things AHF
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

OBJS = ahf_io.o \
       ahf_halos_sfc.o \
       ahf_halos.o \
       ahf_io_sql.o\
       patchtree2halos.o

# Catch empty make call
all:
	${MAKE} libahf2.a

# Cleaning up
clean:
	rm -f libahf2.a $(OBJS)

# The library itself
libahf2.a: $(OBJS)
	${AR} -r libahf2.a $(OBJS)

# The individual routines
ahf_halos.o: ahf_halos.c
	$(CC) $(CFLAGS) -c -o ahf_halos.o ahf_halos.c

ahf_halos_sfc.o: ahf_halos_sfc.c ahf_halos_sfc.h
	$(CC) $(CFLAGS) -c -o ahf_halos_sfc.o ahf_halos_sfc.c

ahf_io.o: ahf_io.c ahf_io.h
	$(CC) $(CFLAGS) -c -o ahf_io.o ahf_io.c

ahf_io_sql.o: ahf_io_sql.c ahf_io_sql.h ahf_io_sql_adt.h
	$(CC) $(CFLAGS) -c -o ahf_io_sql.o ahf_io_sql.c

patchtree2halos.o: patchtree2halos.c
	$(CC) $(CFLAGS) -c patchtree2halos.c
