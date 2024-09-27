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

.PHONY: clean veryclean all

#*******************************************************************#
DEFINEFLAGS = ${MASTER_DEFINEFLAGS} -DAHF

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
	libamr_serial/libamr_serial.a \
	libgravity/libgravity.a \
	libio/libio.a \
	libsfc/libsfc.a \
	libio_serial/libio_serial.a \
	libutility/libutility.a

LIBAHF  = libahf/libahf.a

SUBDIRS = libahf libgravity libio_serial libutility libamr_serial libio libsfc
.PHONY: $(SUBDIRS)

all: AHF


# what's to make ?
# ----------------
AHF:	main.o ${OBJS} $(SUBDIRS)
	$(CC) $(CFLAGS) ${PROFILE} -o $(BIN)/AHF-$(VERSION)  main.o ${OBJS} ${LIBAHF}  ${LIBS} -l${LIB_MATH} ${LFLAGS}


clean:
	rm -f *.o AMIGA* AHF* CONVERT *.il
	$(MAKE) -C libsfc        clean
	$(MAKE) -C libio         clean
	$(MAKE) -C libio_serial  clean
	$(MAKE) -C libahf        clean
	$(MAKE) -C libgravity    clean
	$(MAKE) -C libamr_serial clean
	$(MAKE) -C libutility    clean

veryclean:
	rm -f *.o AHF* *~ *~.* *.il *.dSYM
	$(MAKE) -C libsfc        clean
	$(MAKE) -C libio         clean
	$(MAKE) -C libio_serial  clean
	$(MAKE) -C libahf        clean
	$(MAKE) -C libgravity    clean
	$(MAKE) -C libamr_serial clean
	$(MAKE) -C libutility    clean


# dependencies for include files
# ------------------------------
common.h:	param.h\
				tdef.h
tdef.h: 		param.h



# how to make the libraries?
# ------------------------------

$(SUBDIRS):
	$(MAKE) -C $@


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
