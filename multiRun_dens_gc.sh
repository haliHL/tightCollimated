#This bash file is to run multiple simulations with respect to 
#--different average flux
#--different gammac
#for the tight collimated beam laser
#using cNumber theory.

iFile=input.txt
controlType="dens_1000_gc_0.002_0.002_20_randomPhase"
nMax1=1
nMax2=20
init_dens=1000
interval_dens=100
init_gc=0.002
interval_gc=0.002

for ((i=0; i<nMax1; i+=1)) do 
for ((j=0; j<nMax2; j+=1)) do 

dens=$(echo "$init_dens + $interval_dens * $i" | bc -l)
gc=$(echo "$init_gc + $interval_gc * $j" | bc -l | awk '{printf "%.3f", $0}')

printf "tMax 100
nStore 1000
nTrajectory 200
nBin 20
density $dens
gc $gc
controlType $controlType
name dens${dens}_gc${gc}
fast 1" > $iFile

./tightCollimated -f $iFile

done

number=$((1+$i))
echo "Run ${number} of" $nMax1

done

cp $iFile $controlType/"inputControl.txt"