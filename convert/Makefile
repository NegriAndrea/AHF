.PHONY: clean veryclean

simu2tipsy:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
ahf2binary:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
AHFbinary2ascii:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
AHFhalos2geom:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
ramses2gadget:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}

CFLAGS  =       ${OPTIMIZE} ${CCFLAGS} ${DEFINEFLAGS}
LFLAGS  =       ${LNFLAGS}

export DEFINEFLAGS

LIBS = 	../src/libio_serial/libio_serial.a \
	../src/libahf/libahf.a \
	../src/libgravity/libgravity.a \
	../src/libutility/libutility.a \
	../src/libamr_serial/libamr_serial.a


# HOW TO MAKE THE ACTUAL PROGRAMMES
simu2tipsy:		simu2tipsy.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o simu2tipsy simu2tipsy.c ${LIBS} ${LIBS} -lm

AHFhalos2geom:		halos2geom.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o halos2geom halos2geom.c ${LIBS} ${LIBS} -lm

ahf2binary:		ahf2binary.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o ahf2binary ahf2binary.c ${LIBS} ${LIBS} -lm

AHFbinary2ascii:		AHFbinary2ascii.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o AHFbinary2ascii AHFbinary2ascii.c ${LIBS} ${LIBS} -lm

ramses2gadget:		ramses2gadget.f90
	$(FC) $(CFLAGS) $(LFLAGS) -x f95-cpp-input -o ramses2gadget ramses2gadget.f90


# HOW TO MAKE THE LIBRARIES
../src/libio_serial/libio_serial.a:
	cd ../src/libio_serial ; ${MAKE} libio_serial.a

../src/libahf/libahf.a:
	cd ../src/libahf ; ${MAKE} libahf.a

../src/libgravity/libgravity.a:
	cd ../src/libgravity ; ${MAKE} libgravity.a

../src/libamr_serial/libamr_serial.a:
	cd ../src/libamr_serial ; ${MAKE} libamr_serial.a

../src/libutility/libutility.a:
	cd ../src/libutility ; ${MAKE} libutility.a

clean:
	rm -f *.o

veryclean:
	rm -f *.o *~ *~.*
