.PHONY: clean veryclean

ahfCalcHaloID:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
ahfCheckHaloIDs:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
ahfCalcMinStar:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
ahfExtractHalosPdens:	DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
ahfExtractMhalo:	DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
ahfExtractProfiles:	DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
ahfRemUnbound:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
ahfSortHalos:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
simuExtractHalos:	DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
HaloShape:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
cosmology:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
virial:			DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
age:			DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
pmass:			DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}

CFLAGS  =       ${OPTIMIZE} ${CCFLAGS} ${DEFINEFLAGS}
LFLAGS  =       ${LNFLAGS}

export DEFINEFLAGS

LIBS = 	../src/libio_serial/libio_serial.a \
	../src/libahf/libahf.a \
	../src/libgravity/libgravity.a \
	../src/libutility/libutility.a \
	../src/libamr_serial/libamr_serial.a


# HOW TO MAKE THE ACTUAL PROGRAMMES
simuExtractHalos:	simuExtractHalos.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/simuExtractHalos simuExtractHalos.c ${LIBS} ${LIBS} -lm

ahfExtractHalosPdens:	ahfExtractHalosPdens ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/ahfExtractHalosPdens ahfExtractHalosPdens.c ${LIBS} ${LIBS} -lm

virial:	virial.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/virial virial.c ${LIBS} ${LIBS} -lm

age:	age.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/age age.c ${LIBS} ${LIBS} -lm

agelist:	agelist.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/agelist agelist.c ${LIBS} ${LIBS} -lm

cosmology:	cosmology.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/cosmology cosmology.c ${LIBS} ${LIBS} -lm

pmass:		pmass.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/pmass pmass.c -lm

ahfExtractMhalo:		ahfExtractMhalo.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/ahfExtractMhalo ahfExtractMhalo.c -lm

RemoveUnbound:		RemoveUnbound.c
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/RemoveUnbound RemoveUnbound.c -lm

ahfExtractProfiles:		ahfExtractProfiles.c ${LIBS}
		$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/ahfExtractProfiles ahfExtractProfiles.c -lm

HaloShape:		HaloShape.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/HaloShape HaloShape.c ${LIBS} ${LIBS} -lm

ahfCalcMinStar:		ahfCalcMinStar.c
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/ahfCalcMinStar ahfCalcMinStar.c -lm

ahfSortHalos:		ahfSortHalos.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/ahfSortHalos ahfSortHalos.c ${LIBS} ${LIBS} -lm

ahfCheckHaloIDs:		ahfCheckHaloIDs.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/ahfCheckHaloIDs ahfCheckHaloIDs.c ${LIBS} ${LIBS} -lm

ahfCalcHaloID:		ahfCalcHaloID.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/ahfCalcHaloID ahfCalcHaloID.c ${LIBS} ${LIBS} -lm


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
	rm -rf  *.dSYM
	rm -f *.o *~ *~.*
