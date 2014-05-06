#******************************************************************#
#                                                                  #
# Makefile for               AHF                                   #
#                           =====                                  #
#                                                                  #
# version 1.0 (30/01/2005), Alexander Knebe                        #
#                                                                  #
# AHF (AMIGA-Halo-Finder) added by Stuart Gill (30/01/2005)        #
#                                                                  #
#******************************************************************#

.PHONY: clean veryclean

#*******************************************************************#
AHF:			DEFINEFLAGS = ${MASTER_DEFINEFLAGS} -DAHF
AHF2:			DEFINEFLAGS = ${MASTER_DEFINEFLAGS} -DAHF2

export DEFINEFLAGS

# libraries in use
# ----------------
LIB_MATH 	=	m

# flags for compiling all individual routines
# -------------------------------------------
CFLAGS	=	${OPTIMIZE} ${CCFLAGS} ${DEFINEFLAGS}

# flags for linking all parts together
# ------------------------------------
LFLAGS	=	${LNFLAGS}

# all object files
# ----------------
OBJS	=	\
	comm.o \
	common.o\
	startrun.o
  
LIBS =	\
			libgravity/libgravity.a \
			libio/libio.a \
			libsfc/libsfc.a \
			libtree/libtree.a \
			libamr_serial/libamr_serial.a \
			libio_serial/libio_serial.a \
			libutility/libutility.a

LIBAHF  = libahf/libahf.a
LIBAHF2 = libahf2/libahf2.a


# what's to make ?
# ----------------
AHF:	main.o ${OBJS} ${LIBS} ${LIBAHF}
	$(CC) $(CFLAGS) ${PROFILE} -o AHF  main.o ${OBJS} ${LIBAHF}  ${LIBS} -l${LIB_MATH} ${LFLAGS}

AHF2:	main.o ${OBJS} ${LIBS} ${LIBAHF2}
	$(CC) $(CFLAGS) ${PROFILE} -o AHF2 main.o ${OBJS} ${LIBAHF2} ${LIBS} -l${LIB_MATH} ${LFLAGS}


clean:
	rm -f *.o AMIGA* AHF* CONVERT *.il
	cd libsfc ; ${MAKE} clean
	cd libio ; ${MAKE} clean
	cd libio_serial ; ${MAKE} clean
	cd libahf ; ${MAKE} clean
	cd libahf2 ; ${MAKE} clean
	cd libgravity ; ${MAKE} clean
	cd libtree ; ${MAKE} clean
	cd libamr_serial ; ${MAKE} clean
	cd libutility ; ${MAKE} clean

veryclean:
	rm -f *.o AHF* *~ *~.* *.il
	cd libsfc ; ${MAKE} clean
	cd libio ; ${MAKE} clean
	cd libio_serial ; ${MAKE} clean
	cd libahf ; ${MAKE} clean
	cd libahf2 ; ${MAKE} clean
	cd libgravity ; ${MAKE} clean
	cd libtree ; ${MAKE} clean
	cd libamr_serial ; ${MAKE} clean
	cd libutility ; ${MAKE} clean


# dependencies for include files
# ------------------------------
common.h:	param.h\
				tdef.h
tdef.h: 		param.h



# how to make the libraries?
# ------------------------------
libsfc/libsfc.a:
	cd libsfc ; ${MAKE} libsfc.a

libio/libio.a:
	cd libio ; ${MAKE} libio.a

libio_serial/libio_serial.a:
	cd libio_serial ; ${MAKE} libio_serial.a

libahf/libahf.a:
	cd libahf ; ${MAKE} libahf.a

libahf2/libahf2.a:
	cd libahf2 ; ${MAKE} libahf2.a

libgravity/libgravity.a:
	cd libgravity ; ${MAKE} libgravity.a

libtree/libtree.a:
	cd libtree ; ${MAKE} libtree.a

libamr_serial/libamr_serial.a:
	cd libamr_serial ; ${MAKE} libamr_serial.a

libutility/libutility.a:
	cd libutility ; ${MAKE} libutility.a



#dependencies for .o files
# ------------------------
comm.o:	comm.c
	${CC} ${CFLAGS} -c comm.c

common.o:	common.c
	${CC} ${CFLAGS} -c common.c

main.o:			main.c
	${CC} ${CFLAGS} -c main.c

startrun.o:		startrun.c
	${CC} ${CFLAGS} -c startrun.c
