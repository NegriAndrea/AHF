#!/bin/csh

set DATA = /data4/gustavo/MUSIC_PLANCK/DMONLY/HIGH_RES/

foreach GALAXY (`\ls -d $DATA/NewMD*`)
    echo $GALAXY

    set OUTDIR = `basename $GALAXY`
    if(-e $OUTDIR)then
        echo $OUTDIR does exit
    else
        echo $OUTDIR does not exist: just created it...
        mkdir $OUTDIR
    endif

    foreach INFILE (`\ls $GALAXY/*_halos`)
        echo $INFILE
	
        set snapid = `echo $INFILE | awk -F_200c '{print $1}'`
        set snapid = `printf $snapid | tail -c 3`

	./unique_haloids_halos $INFILE $snapid

	mv -f *_halos_uniquehaloids $OUTDIR
    end
end
