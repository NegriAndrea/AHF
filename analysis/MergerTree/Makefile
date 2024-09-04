##################################################################################
# simple (top-level) Makefile for all things MergerTree related
##################################################################################

MergerTree:
	cd src/; \
	make MergerTree; \
	mv -f MergerTree ../bin

MergerRates:
	cd analysis/; \
	make MergerRates; \
	mv -f MergerRates ../bin

clean:
	cd src/; \
	make clean; \
	cd ../analysis/; \
	make clean

veryclean:
	cd src/; \
	make clean; \
	cd ../analysis/; \
	make clean; \
	cd ../; \
	rm -f bin/*

