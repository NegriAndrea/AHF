#!/bin/csh

foreach GALAXY (`\ls -d ./g1.05e11*`)
    foreach INFILE (`\ls $GALAXY/*_halos`)
        echo $INFILE
        set OUTFILE = $INFILE\_uniquehaloids
        
        rm -f bla1 bla2
        touch bla1 bla2
        
        set snapid = `echo $INFILE | awk -F.z '{print $1}'`
        set snapid = `printf $snapid | tail -c 5`

        set HEADER = `head -1 $INFILE`
        echo $HEADER > $OUTFILE
        set NHALOS = `wc -l $INFILE | awk '{print $1}'`
        @ NHALOS = $NHALOS + 1
        set i = 2
        while ($i < $NHALOS)
            set line = `head -$i $INFILE | tail -1`
            
            set haloid       = `echo $line | awk '{print $1}'`
            set hosthaloid   = `echo $line | awk '{print $2}'`
                        
            set uniquehaloid = `expr "$snapid" \* 1000000000000 + "$haloid" + 1`
            if($hosthaloid != -1)then
                set uniquehosthaloid = `expr "$snapid" \* 1000000000000 + "$hosthaloid" + 1`
            else
                set uniquehosthaloid = 0
            endif
            
            echo $uniquehaloid $uniquehosthaloid    >> bla1
            echo $line | awk '{$1=$2=""; print $0}' >> bla2
                        
            @ i = $i + 1
        end
    
        paste bla1 bla2 >> $OUTFILE
        #mv -f $OUTFILE $INFILE
    end
end

rm -f bla1 bla2
