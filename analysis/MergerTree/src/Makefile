#CC		=	gcc
#FLAGS	=	-O2 -fopenmp
LIBS	=	-lm
OBJECTS	=	libio.o libmerit.o libmtree.o libpidmap.o libsnapskipping.o libutil.o

MergerTree: main.c common.h define.h tdef.h $(OBJECTS)
	$(CC) $(CCFLAGS) -o $(BIN)/MergerTree main.c $(OBJECTS) $(LIBS) $(LNFLAGS)


libio.o: libio.c
	$(CC) $(CCFLAGS) -c libio.c $(LIBS)

libmerit.o: libmerit.c
	$(CC) $(CCFLAGS) -c libmerit.c $(LIBS)

libmtree.o: libmtree.c
	$(CC) $(CCFLAGS) -c libmtree.c $(LIBS)

libpidmap.o: libpidmap.c
	$(CC) $(CCFLAGS) -c libpidmap.c $(LIBS)

libsnapskipping.o: libsnapskipping.c
	$(CC) $(CCFLAGS) -c libsnapskipping.c $(LIBS)

libutil.o: libutil.c libutil.h
	$(CC) $(CCFLAGS) -c libutil.c $(LIBS)


.PHONY: clean veryclean

clean:
	rm -f *.o

veryclean:
	rm -f *.o
