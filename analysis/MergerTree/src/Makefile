CC		=	gcc
FLAGS	=	-O2 -fopenmp
LIBS	=	-lm
OBJECTS	=	libio.o libmerit.o libmtree.o libpidmap.o libsnapskipping.o libutil.o

MergerTree: main.c common.h define.h tdef.h $(OBJECTS)
	$(CC) $(FLAGS) -o MergerTree main.c $(OBJECTS) $(LIBS)


libio.o:							libio.c
	$(CC) $(FLAGS) -c libio.c $(LIBS)

libmerit.o:        		libmerit.c
	$(CC) $(FLAGS) -c libmerit.c $(LIBS)

libmtree.o:        		libmtree.c
	$(CC) $(FLAGS) -c libmtree.c $(LIBS)

libpidmap.o:        	libpidmap.c
	$(CC) $(FLAGS) -c libpidmap.c $(LIBS)

libsnapskipping.o:    libsnapskipping.c
	$(CC) $(FLAGS) -c libsnapskipping.c $(LIBS)

libutil.o:      			libutil.c libutil.h
	$(CC) $(FLAGS) -c libutil.c $(LIBS)



clean:
	rm -f *.o
