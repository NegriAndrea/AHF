.PHONY: clean veryclean

DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
CFLAGS  =       ${OPTIMIZE} ${CCFLAGS} ${DEFINEFLAGS}
LFLAGS  =       ${LNFLAGS}

export DEFINEFLAGS

LIBS = 	../src/libio_serial/libio_serial.a \
	../src/libahf/libahf.a \
	../src/libgravity/libgravity.a \
	../src/libutility/libutility.a \
	../src/libamr_serial/libamr_serial.a

DIRLIBS = 	../src/libio_serial \
	../src/libahf \
	../src/libgravity \
	../src/libutility \
	../src/libamr_serial
.PHONY: $(DIRLIBS)

# HOW TO MAKE THE ACTUAL PROGRAMMES
simu2tipsy:		simu2tipsy.c ${DIRLIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o simu2tipsy simu2tipsy.c ${LIBS} -lm

AHFhalos2geom:		AHFhalos2geom.c ${DIRLIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o halos2geom AHFhalos2geom.c ${LIBS} -lm

ahf2binary:		ahf2binary.c ${DIRLIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o ahf2binary ahf2binary.c ${LIBS} -lm

AHFbinary2ascii:		AHFbinary2ascii.c ${DIRLIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o AHFbinary2ascii AHFbinary2ascii.c ${LIBS} -lm

ramses2gadget:		ramses2gadget.f90
	$(FC) $(CFLAGS) $(LFLAGS) -x f95-cpp-input -o ramses2gadget ramses2gadget.f90


## HOW TO MAKE THE LIBRARIES
$(DIRLIBS):
	$(MAKE) -C $@

clean:
	rm -f *.o

veryclean:
	rm -f *.o *~ *~.*
