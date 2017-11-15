##################################################################################
# WHATEVER YOU DO, PLEASE DO IT TO Makefile.config AND LEAVE THIS MAKEFILE ALONE #
##################################################################################
include Makefile.config


#------------------------------------------------------------------#
# the menu                                                         #
#------------------------------------------------------------------#
MASTER_DEFINEFLAGS	=	$(DEFINEFLAGS)

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

# everything in src/
#===================
AHF:	FORCE dirs
	cd src;\
	${MAKE} AHF;\
	mv -f AHF ../bin/AHF-v1.0-093

AHF2:	FORCE dirs
	cd src;\
	${MAKE} AHF2;\
	mv -f AHF2 ../bin/AHF-v2.0-000


# everything in convert/
#=======================
simu2tipsy:	FORCE dirs
		cd convert;\
		${MAKE} simu2tipsy;\
		mv -f simu2tipsy ../bin

ahf2binary:	FORCE dirs
		cd convert;\
		${MAKE} ahf2binary;\
		mv -f ahf2binary ../bin

AHFbinary2ascii:	FORCE dirs
		cd convert;\
		${MAKE} AHFbinary2ascii;\
		mv -f AHFbinary2ascii ../bin

AHFhalos2geom:	FORCE dirs
		cd convert;\
		${MAKE} halos2geom;\
		mv -f halos2geom ../bin

ramses2gadget:	FORCE dirs
		cd convert;\
		${MAKE} ramses2gadget;\
		mv -f ramses2gadget ../bin

Hdf2GADGET:	FORCE 
		echo "";echo "please change to convert/Hdf2GADGET/ and refer to the README in there...";echo ""

Hdf2GADGET:	FORCE 
		echo "";echo "please change to convert/tostd/ and refer to the README in there...";echo ""

# everything in tools/
#=====================
ahfCalcHaloID:	FORCE dirs
		cd tools;\
		${MAKE} ahfCalcHaloID;\
		mv -f ahfCalcHaloID ../bin

ahfCheckHaloIDs:	FORCE dirs
		cd tools;\
		${MAKE} ahfCheckHaloIDs;\
		mv -f ahfCheckHaloIDs ../bin

ahfSortHalos:	FORCE dirs
		cd tools;\
		${MAKE} ahfSortHalos;\
		mv -f ahfSortHalos ../bin

simuExtractHalos:	FORCE dirs
		cd tools;\
		${MAKE} simuExtractHalos;\
		mv -f simuExtractHalos ../bin

ahfExtractHalosPdens:	FORCE dirs
		cd tools;\
		${MAKE} ahfExtractHalosPdens;\
		mv -f ahfExtractHalosPdens ../bin

virial:	FORCE dirs
		cd tools;\
		${MAKE} virial;\
		mv -f virial ../bin

age:	FORCE dirs
		cd tools;\
		${MAKE} age;\
		mv -f age ../bin

agelist:	FORCE dirs
		cd tools;\
		${MAKE} agelist;\
		mv -f agelist ../bin

cosmology:	FORCE dirs
		cd tools;\
		${MAKE} cosmology;\
		mv -f cosmology ../bin

pmass:		FORCE dirs
		cd tools;\
		${MAKE} pmass;\
		mv -f pmass ../bin

RemoveUnbound:	FORCE dirs
		cd tools;\
		${MAKE} RemoveUnbound;\
		mv -f RemoveUnbound ../bin

ahfExtractMhalo:		FORCE dirs
		cd tools;\
		${MAKE} ahfExtractMhalo;\
		mv -f ahfExtractMhalo ../bin

ahfExtractProfiles:		FORCE dirs
		cd tools;\
		${MAKE} ahfExtractProfiles;\
		mv -f ahfExtractProfiles ../bin

HaloShape:		FORCE dirs
		cd tools;\
		${MAKE} HaloShape;\
		mv -f HaloShape ../bin

ahfCalcMinStar:		FORCE dirs
		cd tools;\
		${MAKE} ahfCalcMinStar;\
		mv -f ahfCalcMinStar ../bin



# everything in analysis/
#========================
simuPdens:	FORCE dirs
		cd analysis;\
		${MAKE} simuPdens;\
		mv -f simuPdens ../bin

simuCroCo:	FORCE dirs
		cd analysis;\
		${MAKE} simuCroCo;\
		mv -f simuCroCo ../bin

simuPk:	FORCE dirs
		cd analysis;\
		${MAKE} simuPk;\
		mv -f simuPk ../bin

simuXi:	FORCE dirs
		cd analysis;\
		${MAKE} simuXi;\
		mv -f simuXi ../bin

simuVstat:	FORCE dirs
		cd analysis;\
		${MAKE} simuVstat;\
		mv -f simuVstat ../bin

simuSigmaR:	FORCE dirs
		cd analysis;\
		${MAKE} simuSigmaR;\
		mv -f simuSigmaR ../bin

MergerRates:			FORCE dirs
		cd analysis;\
		${MAKE} MergerRates;\
		mv -f MergerRates ../bin

MergerTree:			FORCE dirs
		cd analysis;\
		${MAKE} MergerTree;\
		mv -f MergerTree ../bin

MergerTreeMPI:			FORCE dirs
		cd analysis;\
		${MAKE} MergerTreeMPI;\
		mv -f MergerTreeMPI ../bin

ahfFindHalo:			FORCE dirs
		cd analysis;\
		${MAKE} ahfFindHalo;\
		mv -f ahfFindHalo ../bin

ahfFindHaloPairs:			FORCE dirs
		cd analysis;\
		${MAKE} ahfFindHaloPairs;\
		mv -f ahfFindHaloPairs ../bin

ahfFindCrossPairs:			FORCE dirs
		cd analysis;\
		${MAKE} ahfFindCrossPairs;\
		mv -f ahfFindCrossPairs ../bin

ahfHaloHistory:			FORCE dirs
		cd analysis;\
		${MAKE} ahfHaloHistory;\
		mv -f ahfHaloHistory ../bin

ahfSubhaloAccretion:			FORCE dirs
		cd analysis;\
		${MAKE} ahfSubhaloAccretion;\
		mv -f ahfSubhaloAccretion ../bin

ahfSubhaloAccretionStats:			FORCE dirs
		cd analysis;\
		${MAKE} ahfSubhaloAccretionStats;\
		mv -f ahfSubhaloAccretionStats ../bin

ahfSubCheck:		FORCE dirs
		cd analysis;\
		${MAKE} ahfSubCheck;\
		mv -f ahfSubCheck ../bin

ahfXi:			FORCE dirs
		cd analysis;\
		${MAKE} ahfXi;\
		mv -f ahfXi ../bin

sigmaH:			FORCE dirs
		cd analysis;\
		${MAKE} sigmaH;\
		mv -f sigmaH ../bin


#-------------------------------------------------------------------#
# "make clean" 
#-------------------------------------------------------------------#
clean:	FORCE
	@echo '';\
	echo '*=======================================================================*';\
	echo ' ${MAKE} clean: removing all unneeded files';\
	echo '*=======================================================================*';\
	echo '';\
	echo ''
	cd src; ${MAKE} clean;\
	cd ../convert; ${MAKE} clean;\
	cd ../tools; ${MAKE} clean;\
	cd ../analysis; ${MAKE} clean

#-------------------------------------------------------------------#
# "make veryclean" 
#-------------------------------------------------------------------#
veryclean:	FORCE
	@echo '';\
	echo '*=======================================================================*';\
	echo ' ${MAKE} veryclean: removing all unneeded files (incl. emacs backup files!)';\
	echo '*=======================================================================*';\
	echo '';\
	echo ''
	rm -f *~ *~.* bin/*;\
	cd src; ${MAKE} veryclean;\
	cd ../convert; ${MAKE} veryclean;\
	cd ../tools; ${MAKE} veryclean;\
	cd ../analysis; ${MAKE} veryclean

FORCE:	

dirs:
	mkdir -p bin
