CC	=	gcc
FLAGS	=	-O2 -fopenmp
LIBS	=	-lm
OBJECTS	=	libutil.o

MergerRates: MergerRates.c $(OBJECTS)
	$(CC) $(FLAGS) -o MergerRates MergerRates.c $(OBJECTS) $(LIBRARIES)

libutil.o:	../src/libutil.c ../src/libutil.h
	$(CC) $(FLAGS) -c libutil.o ../src/libutil.c $(LIBS)

clean:
	rm -f *.o 
