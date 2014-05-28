# If ${OPTIMIZE}, ${CCFLAGS} and/or ${DEFINEFLAGS} are set, they will be
# used for the compilation
CFLAGS+=${OPTIMIZE} ${CCFLAGS} ${DEFINEFLAGS}
CFLAGS+=-Wall -Winline -fopenmp
# If a different MAKE is specified use that one, otherwise try your luck
# with just make; maybe it will work.
MAKE?=make

# The same for ar
AR?=ar

# PHONY targets
.PHONY: all clean

## Modified by F.Campos: subcube.o move to the first position trying to avoid 
## linking error: undefined reference to `table_find_subcube' in generate_tree function
OBJS = generate_tree.o cubekey.o subcube.o patch.o utilities.o set_patch_pos.o set_patch_radii.o set_patch_trunk.o

# Catch empty make call
all:
	${MAKE} libtree.a

# Cleaning up
clean:
	rm -f libtree.a $(OBJS) *.h.gch ck2pos.o ck2pos

# The library itself
libtree.a: $(OBJS)
	${AR} -r libtree.a $(OBJS)

# The individual routines stuff
subcube.o:	subcube.c
	$(CC) $(CFLAGS) -c subcube.c

generate_tree.o:	generate_tree.c
	$(CC) $(CFLAGS) -c generate_tree.c 

cubekey.o:	cubekey.c
	$(CC) $(CFLAGS) -c cubekey.c

patch.o:	patch.c
	$(CC) $(CFLAGS) -c patch.c

set_patch_pos.o:	set_patch_pos.c
	$(CC) $(CFLAGS) -c set_patch_pos.c

set_patch_radii.o:	set_patch_radii.c
	$(CC) $(CFLAGS) -c set_patch_radii.c

set_patch_trunk.o:	set_patch_trunk.c
	$(CC) $(CFLAGS) -c set_patch_trunk.c

utilities.o:	utilities.c utilities.h
	$(CC) $(CFLAGS) -c utilities.c

