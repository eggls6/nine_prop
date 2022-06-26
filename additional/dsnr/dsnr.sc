#!/bin/bash
# this script calcualtes Deflection Singal to Noise Ratios from Radau_sml.kin files
# needs gnuplot > v4.2 

outfile="dsnr.txt"

rm $outfile

file="Radau_sml.kin"


nscen=5
ncase=3 #randkick, surekick, nokick
tname="Itokawa"

dv='($9**2+$10**2+$11**2) prefix dv'
vast='($6**2+$7**2+$8**2) prefix vast'
dvsure='($9**2+$10**2+$11**2) prefix dvsure'
vastnom='($6**2+$7**2+$8**2) prefix vastnom'


casename=("randkick" "surekick" "nokick")

for (( i=0; i<$nscen; i++ )); do # different dates

 cd "scenario$i"

 cd $tname
 
 #for (( j=0; j<$ncase; j++ )); do #rnd sure no
 
 echo "# random kick
 #get sigma v of asteroid at kicktime
 stat '${casename[$1]}/$file' $vast
 #get stats on the deflection change in the asteroids velocity
 stat '${casename[$1]}/$file' $dv

 #sure kick
 #get sigma v of asteroid at kicktime
 stat '${casename[$2]}/$file' $vastnom
 #get stats on the deflection change in the asteroids velocity
 stat '${casename[$2]}/$file' $dvsure
 
 
 dsnrnom=dvsure_median/vastnom_stddev
 dsnrmax=dv_max/vast_stddev
 dsnrmin=dv_min/vast_stddev
 
 sprintf('scenario $i:     %g   %g  %g',dsnrmax,dsnrmin,dsnrnom)    
 " > dsnr.gnu
 gnuplot dsnr.gnu >> ../../../$outfile
 
 done
done

exit 0
 
 
 
 