#This bash file is to run multiple simulations with respect to 
#--different average flux
#--different gammac
#for the beam laser
#using cNumber theory.

iFile=input.txt
controlType="tau1.0_dens_100_100_10_gc_0.005_0.005_10_pois0"
nMax1=10
nMax2=10
init_dens=100
interval_dens=100
init_gc=0.005
interval_gc=0.005

for ((i=0; i<nMax1; i+=1)) do 
for ((j=0; j<nMax2; j+=1)) do 

dens=$(echo "$init_dens + $interval_dens * $i" | bc -l)
gc=$(echo "$init_gc + $interval_gc * $j" | bc -l | awk '{printf "%.4f", $0}')
kappa=$(echo "100.0 / $gc" |bc -l)

printf "dt 0.01
tMax 100
nStore 1000
nTrajectory 200
nBin 20
yWall 1.0
lambda 1.0
deltaZ 0
deltaPz 0
tau 1.0
density $dens
mAtom 1
rabi 10
kappa $kappa
detuning 0
T2 0
controlType $controlType
name tau1.0_dens${dens}_gc${gc}
pois 0
fast 1" > $iFile

./beamLaser -f $iFile

done

number=$((1+$i))
echo "Run ${number} of" $nMax1

done

cp $iFile $controlType/"inputControl.txt"