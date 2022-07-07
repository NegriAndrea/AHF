#!/bin/csh

foreach GALAXY (`\ls -d ./g*`)
    foreach INFILE (`\ls $GALAXY/*_particles`)
        echo $INFILE
        set snapid = `echo $INFILE | awk -F.z '{print $1}'`
        set snapid = `printf $snapid | tail -c 5`
        ./unique_haloids_particles $INFILE $snapid
    end
end
