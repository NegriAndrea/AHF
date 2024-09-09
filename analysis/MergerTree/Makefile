##################################################################################
# simple (top-level) Makefile for all things MergerTree related
##################################################################################

.PHONY: MergerTree MergerRates clean veryclean

MergerTree:
	make -C src MergerTree

MergerRates:
	make -C analysis MergerRates

clean:
	make -C src clean
	make -C analysis clean

veryclean:
	make -C src veryclean
	make -C analysis veryclean
	rm -f bin/*

