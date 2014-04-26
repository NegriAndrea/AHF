#!/bin/csh

set FILE      = "lcdm/snapshot_057.z0.000.AHF_halos"
set OUTFILE   = `basename $FILE`
set OUTFILE   = "$OUTFILE.AHFascii"

set HEADER    = "$FILE"
set MULTIMASS = "0"
set BOXSIZE   = "250"
set NPART     = `wc -l $FILE | awk '{print $1}'`
set NPART     = `expr $NPART`
@ NPART       = $NPART - 1
set NSPECIES  = "1"
set TOTALMASS = "1e20"
set NSTEPS    = "1000"
set OMEGA0    = "0.27"
set LAMBDA0   = "0.73"
set PMASS     = "1"
set ZINIT     = "200"
set ZNOW      = `echo $FILE | awk -F.z '{print $2}' | awk -F.AHF '{print $1}'`

echo "# Header = "           $HEADER     > $OUTFILE
echo "# MultiMass = "        $MULTIMASS >> $OUTFILE
echo "# Boxsize = "          $BOXSIZE   >> $OUTFILE
echo "# NumberOfPart = "     $NPART     >> $OUTFILE
echo "# NumberOfSpecies = "  $NSPECIES  >> $OUTFILE
echo "# TotalMass = "        $TOTALMASS >> $OUTFILE
echo "# NumberOfTimestep = " $NSTEPS    >> $OUTFILE
echo "# Omega0 = "           $OMEGA0    >> $OUTFILE
echo "# Lambda0 = "          $LAMBDA0   >> $OUTFILE
echo "# PMass = "            $PMASS     >> $OUTFILE
echo "# ZInitial = "         $ZINIT     >> $OUTFILE
echo "# ZCurrent = "         $ZNOW      >> $OUTFILE

sed '1, 1d' $FILE | awk '{print $6/1000.,$7/1000.,$8/1000.,$10,$11,$12,$4}' >> $OUTFILE
