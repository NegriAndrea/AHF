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

OBJS = 	coarse_to_fine.o \
	fine_to_coarse.o\
	generate_grids.o\
	get_nnodes.o\
	go_down.o\
	refine_grid.o\
	density.o \
	lltools.o\
	relink.o

# Catch empty make call
all:
	${MAKE} libamr_serial.a

# Cleaning up
clean:
	rm -f libamr_serial.a $(OBJS)

# The library itself
libamr_serial.a: $(OBJS)
	${AR} -r libamr_serial.a $(OBJS)


# The individual routines stuff
coarse_to_fine.o: coarse_to_fine.c
	$(CC) $(CFLAGS) -c -o coarse_to_fine.o coarse_to_fine.c

fine_to_coarse.o: fine_to_coarse.c
	$(CC) $(CFLAGS) -c -o fine_to_coarse.o fine_to_coarse.c

generate_grids.o: generate_grids.c
	$(CC) $(CFLAGS) -c -o generate_grids.o generate_grids.c

get_nnodes.o: get_nnodes.c
	$(CC) $(CFLAGS) -c -o get_nnodes.o get_nnodes.c

go_down.o: go_down.c
	$(CC) $(CFLAGS) -c -o go_down.o go_down.c

refine_grid.o: refine_grid.c
	$(CC) $(CFLAGS) -c -o refine_grid.o refine_grid.c

density.o: density.c
	$(CC) $(CFLAGS) -c -o density.o density.c

lltools.o: lltools.c
	$(CC) $(CFLAGS) -c -o lltools.o lltools.c

relink.o: relink.c
	$(CC) $(CFLAGS) -c -o relink.o relink.c


