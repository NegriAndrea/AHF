#!/bin/csh
date
setenv OMP_NUM_THREADS 16

set MERGERTREECODE = "MergerTree"  # note, this will also be used as the prefix for the SUSSING2013 output file generated at the end of this script

foreach MODEL ("HiSURFS")

        set DATADIR    = "./$MODEL/"

	######################################
	# generate input file for MergerTree #
	######################################
	rm -f prefix
	touch prefix
	# note, we are assuming here that you ran the serial version of AHF, i.e. only one _particles file per snapshot
	foreach file ( `\ls -r $DATADIR/*_particles | awk -F_particles '{print $1}'` )
		echo $file >> prefix
	end
	\ls -r $DATADIR/*_particles | wc -l > nfiles
	\ls -r $DATADIR/*_particles > particles
	cat nfiles particles prefix > MergerTree.input
	rm -f nfiles particles prefix

	#############################
	# eventually run MergerTree #
	#############################
	./$MERGERTREECODE < MergerTree.input


        #####################################################
        # generate mtree file compliant with Sussing Format #
        #####################################################
        set MERGERTREECODE = `basename $MERGERTREECODE`
        set VERSION   = "1"
	set HEADER    = $MERGERTREECODE
        set MTREEDIR  = $DATADIR
        set OUTFILE   = $MERGERTREECODE\_$MODEL.txt

	rm -f $OUTFILE $OUTFILE.tmp
	touch $OUTFILE
	touch $OUTFILE.tmp

	echo $VERSION > header.txt
	echo $HEADER >> header.txt

	# get number of files
	ls -r $MTREEDIR/*_mtree > ls.txt
	wc -l ls.txt | awk '{print $1}' > nfiles.txt
	rm -f ls.txt

	# get first snapshot
	echo 127 > first.txt

	# get last snapshot
	echo 10 > last.txt

	# [TotNsnapshots]    [FirstSnapshot]    [Lastsnapshot]
	paste nfiles.txt first.txt last.txt > files.txt
	rm -f nfiles.txt first.txt last.txt


	# generate tree file simultaneously counting nhalos
	rm -f nhalos.txt
	touch nhalos.txt
	foreach file ( `ls -r $MTREEDIR/*_mtree` )
	    set FIRSTHALOID = `head -2 $file | tail -1 | awk '{print $1}'`
	    set FIRSTHALOID = `expr $FIRSTHALOID`
	    if ($FIRSTHALOID != 0) then
		echo $file $FIRSTHALOID `head -1 $file`
		head -1 $file >> nhalos.txt
		sed '1, 1 d' $file >> $OUTFILE.tmp
	    endif
	end

	awk '{a+=$0}END{print a}' nhalos.txt | tail -1 >> sumnhalos.txt
	rm -f nhalos.txt

	set nterminated = 0

	# update sumnhalos to account for haloes terminated with this script
	set sumnhalos = `head -1 sumnhalos.txt`
	set sumnhalos = `expr $sumnhalos`
	echo $sumnhalos
	@ sumnhalos = $sumnhalos + $nterminated
	echo $sumnhalos > sumnhalos.txt

	# [TotNhalos]    [TotNsnapshots]    [FirstSnapshot]    [Lastsnapshot]
	paste sumnhalos.txt files.txt >> header.txt
	rm -f sumnhalos.txt files.txt

	# cat header and actual tree together
	cat header.txt $OUTFILE.tmp > $OUTFILE
	echo "END" >> $OUTFILE
	rm -f header.txt
	rm -f $OUTFILE.tmp

end # MODEL-loop

echo FINISHED
date
