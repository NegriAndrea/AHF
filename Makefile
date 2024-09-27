##################################################################################
# WHATEVER YOU DO, PLEASE DO IT TO Makefile.config AND LEAVE THIS MAKEFILE ALONE #
##################################################################################
include Makefile.config


#------------------------------------------------------------------#
# the menu                                                         #
#------------------------------------------------------------------#
MASTER_DEFINEFLAGS	=	$(DEFINEFLAGS)

#------------------------------------------------------------------#
# the bin directory                                                #
#------------------------------------------------------------------#
BIN=$(PWD)/bin
VERSION=v1.0-116
#------------------------------------------------------------------#
# make settings available to all other Makefiles                   #
#------------------------------------------------------------------#
export CC
export FC
export OPTIMIZE
export CCFLAGS
export LNFLAGS
export MASTER_DEFINEFLAGS
export MAKE
export BIN
export VERSION

ALL=AHF simu2tipsy ahf2binary AHFhalos2geom ramses2gadget Hdf2GADGET ahfCalcHaloID ahfCheckHaloIDs ahfSortHalos simuExtractHalos ahfExtractHalosPdens virial age agelist cosmology pmass RemoveUnbound ahfExtractMhalo ahfExtractProfiles HaloShape ahfCalcMinStar simuPdens simuPk simuCroCo simuXi simuVstat simuSigmaR MergerTree ahfFindHalo ahfFindHaloPairs ahfFindCrossPairs ahfHaloHistory ahfSubhaloAccretion ahfSubhaloAccretionStats ahfSubCheck ahfXi sigmaH


.PHONY: ${ALL} veryclean clean cleaner

all: ${ALL}

# everything in src/
#===================
AHF:	bin
	${MAKE} -C src AHF
	#cp src/AHF bin/AHF-v1.0-116
	#mv -f src/AHF bin/AHF-v1.0-116
	echo $(BIN)


# everything in convert/
#=======================
simu2tipsy:	 bin
		${MAKE} -C convert simu2tipsy
		mv -f convert/simu2tipsy bin

ahf2binary:	 bin
		${MAKE} -C convert ahf2binary
		mv -f convert/ahf2binary bin

AHFbinary2ascii:	 bin
		${MAKE} -C convert AHFbinary2ascii
		mv -f convert/AHFbinary2ascii bin

AHFhalos2geom:	 bin
		${MAKE} -C convert AHFhalos2geom
		mv -f convert/halos2geom bin

ramses2gadget:	 bin
		${MAKE} -C convert ramses2gadget
		mv -f convert/ramses2gadget bin

Hdf2GADGET:	 
		echo "";echo "please change to convert/Hdf2GADGET/ and refer to the README in there...";echo ""

Hdf2GADGET:	 
		echo "";echo "please change to convert/tostd/ and refer to the README in there...";echo ""

# everything in tools/
#=====================
ahfCalcHaloID:	 bin
		${MAKE} -C tools ahfCalcHaloID

ahfCheckHaloIDs:	 bin
		${MAKE} -C tools  ahfCheckHaloIDs

ahfSortHalos:	 bin
		${MAKE} -C tools  ahfSortHalos

simuExtractHalos:	 bin
		${MAKE}  -C tools simuExtractHalos

ahfExtractHalosPdens:	 bin
		${MAKE}  -C tools ahfExtractHalosPdens

virial:	 bin
		${MAKE}  -C tools virial

age:	 bin
		${MAKE} -C tools  age

agelist:	 bin
		${MAKE} -C tools  agelist

cosmology:	 bin
		${MAKE} -C tools  cosmology

pmass:		 bin
		${MAKE}  -C tools pmass

RemoveUnbound:	 bin
		${MAKE} -C tools RemoveUnbound

ahfExtractMhalo:		 bin
		${MAKE}  -C tools ahfExtractMhalo

ahfExtractProfiles:		 bin
		${MAKE}  -C tools ahfExtractProfiles

HaloShape:		 bin
		${MAKE}  -C tools HaloShape

ahfCalcMinStar:		 bin
		${MAKE}  -C tools ahfCalcMinStar



# everything in analysis/
#========================
simuPdens:	 bin
		${MAKE} -C analysis simuPdens

simuCroCo:	 bin
		${MAKE} -C analysis simuCroCo

simuPk:	 bin
		${MAKE} -C analysis simuPk

simuXi:	 bin
		${MAKE} -C analysis  simuXi

simuVstat:	 bin
		${MAKE} -C analysis  simuVstat

simuSigmaR:	 bin
		${MAKE} -C analysis  simuSigmaR

MergerTree:			 bin
		${MAKE} -C analysis/MergerTree  MergerTree

ahfFindHalo:			 bin
		${MAKE} -C analysis  ahfFindHalo

ahfFindHaloPairs:			 bin
		${MAKE} -C analysis  ahfFindHaloPairs

ahfFindCrossPairs:			 bin
		${MAKE} -C analysis  ahfFindCrossPairs

ahfHaloHistory:			 bin
		${MAKE} -C analysis  ahfHaloHistory

ahfSubhaloAccretion:			 bin
		${MAKE} -C analysis  ahfSubhaloAccretion

ahfSubhaloAccretionStats:			 bin
		${MAKE} -C analysis  ahfSubhaloAccretionStats

ahfSubCheck:		 bin
		${MAKE} -C analysis  ahfSubCheck

ahfXi:			 bin
		${MAKE} -C analysis  ahfXi

sigmaH:			 bin
		${MAKE} -C analysis  sigmaH


#-------------------------------------------------------------------#
# "make clean" 
#-------------------------------------------------------------------#
clean:	
	@echo '';\
	echo '*=======================================================================*';\
	echo ' ${MAKE} clean: removing all unneeded files';\
	echo '*=======================================================================*';\
	echo '';\
	echo ''
	$(MAKE) -C src                  clean
	$(MAKE) -C convert              clean
	$(MAKE) -C tools                clean
	$(MAKE) -C analysis             clean
	$(MAKE) -C analysis/MergerTree/ clean

#-------------------------------------------------------------------#
# "make veryclean" 
#-------------------------------------------------------------------#
veryclean:	
	@echo '';\
	echo '*=======================================================================*';\
	echo ' ${MAKE} veryclean: removing all unneeded files (incl. emacs backup files!)';\
	echo '*=======================================================================*';\
	echo '';\
	echo ''
	rm -rf  bin/*.dSYM
	rm -f *~ *~.* bin/* *.dSYM
	$(MAKE) -C src                  veryclean
	$(MAKE) -C convert              veryclean
	$(MAKE) -C tools                veryclean
	$(MAKE) -C analysis             veryclean
	$(MAKE) -C analysis/MergerTree/ veryclean


cleaner: veryclean


bin:
	mkdir -p $(BIN)
