simuPdens:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
simuCroCo:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
simuPk:			DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
simuXi:			DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
simuVstat:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
simuSigmaR:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
MergerRates:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
MergerTree:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
MergerTreeMPI:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
ahfSubhaloAccretion:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
ahfSubhaloAccretionStats:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
ahfFindHalo:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
ahfFindHaloPairs:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
ahfFindCrossPairs:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
ahfHaloHistory:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
sigmaH:		DEFINEFLAGS	=	-DAHF -DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}

CFLAGS  =       ${OPTIMIZE} ${CCFLAGS} ${DEFINEFLAGS}
LFLAGS  =       ${LNFLAGS}

export DEFINEFLAGS


OBJS	=	\
	comm.o \
  common.o \
	startrun.o

LIBS =	../src/libahf/libahf.a \
			../src/libgravity/libgravity.a \
			../src/libio/libio.a \
			../src/libsfc/libsfc.a \
			../src/libamr_serial/libamr_serial.a \
			../src/libio_serial/libio_serial.a \
			../src/libutility/libutility.a

# HOW TO MAKE THE ACTUAL PROGRAMMES
#-----------------------------------
simuPdens:	simuPdens.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/simuPdens simuPdens.c ${LIBS} ${LIBS} -lm

simuCroCo:	simuCroCo.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/simuCroCo simuCroCo.c ${LIBS} ${LIBS} -lm

simuPk:	simuPk.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/simuPk simuPk.c ${LIBS} ${LIBS} -lm

simuXi:	simuXi.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/simuXi simuXi.c ${LIBS} ${LIBS} -lm

simuVstat:	simuVstat.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/simuVstat simuVstat.c ${LIBS} ${LIBS} -lm

simuSigmaR:	simuSigmaR.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/simuSigmaR simuSigmaR.c ${LIBS} ${LIBS} -lm

MergerRates:	MergerRates.c  ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/MergerRates  ${LIBS} MergerRates.c  ${LIBS} -lm

MergerTree:	MergerTree.c  ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/MergerTree  ${LIBS} MergerTree.c  ${LIBS} -lm

MergerTreeMPI:	MergerTreeMPI.c  ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/MergerTreeMPI  ${LIBS} MergerTreeMPI.c  ${LIBS} -lm

ahfHaloHistory:	ahfHaloHistory.c  ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/ahfHaloHistory  ${LIBS} ahfHaloHistory.c  ${LIBS} -lm

ahfSubhaloAccretion:	ahfSubhaloAccretion.c  ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/ahfSubhaloAccretion  ${LIBS} ahfSubhaloAccretion.c  ${LIBS} -lm

ahfSubhaloAccretionStats:	ahfSubhaloAccretionStats.c  ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/ahfSubhaloAccretionStats  ${LIBS} ahfSubhaloAccretionStats.c  ${LIBS} -lm

ahfFindHalo:	ahfFindHalo.c
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/ahfFindHalo ahfFindHalo.c -lm

ahfFindHaloPairs:	ahfFindHaloPairs.c
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/ahfFindHaloPairs ahfFindHaloPairs.c -lm

ahfFindCrossPairs:	ahfFindCrossPairs.c
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/ahfFindCrossPairs ahfFindCrossPairs.c -lm

ahfHaloHistoryFortran:		ahfHaloHistory.f
	$(FC) -O2 -o $(BIN)/ahfHaloHistory ahfHaloHistory.f

ahfSubCheck:		ahfSubCheck.f
	$(FC) -O2 -o $(BIN)/ahfSubCheck ahfSubCheck.f

ahfXi:		ahfXi.f
	$(FC) -O2 -o $(BIN)/ahfXi ahfXi.f

sigmaH:	sigmaH.o ${OBJS} ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN)/sigmaH sigmaH.c ${LIBS} ${OBJS} ${LIBS} -lm



# HOW TO MAKE THE LIBRARIES
#---------------------------
../src/libio_serial/libio_serial.a:
	cd ../src/libio_serial ; ${MAKE} libio_serial.a

../src/libio/libio.a:
	cd ../src/libio ; ${MAKE} libio.a

../src/libahf/libahf.a:
	cd ../src/libahf ; ${MAKE} libahf.a

../src/libgravity/libgravity.a:
	cd ../src/libgravity ; ${MAKE} libgravity.a

../src/libsfc/libsfc.a:
	cd ../src/libsfc ; ${MAKE} libsfc.a

../src/libamr_serial/libamr_serial.a:
	cd ../src/libamr_serial ; ${MAKE} libamr_serial.a

../src/libutility/libutility.a:
	cd ../src/libutility ; ${MAKE} libutility.a


# dependencies for .o files
# ------------------------
comm.o:	../src/comm.c
	${CC} ${CFLAGS} -c ../src/comm.c

common.o:	../src/common.c
	${CC} ${CFLAGS} -c ../src/common.c

startrun.o:		../src/startrun.c
	${CC} ${CFLAGS} -c ../src/startrun.c

clean:
	rm -f  *.o

veryclean:
	rm -f  *.o *~ *~.*
