source parameters_operational.bash

firstin=$1              # greastest domain .in file as argument

basein=son_in

j=0
for i in $(echo ${son_grids[@]});do
if [ $j == 0 ];then
  fatherin=$firstin
else
  fatherin=${basein}_${j}.in
fi
j=$(($j + 1))
sonin=${basein}_${j}.in
nestson=$i
echo $fatherin
echo $sonin
echo $nestson
cp $nestson $mercatordata                 #moving son grid to mercator directory for interp  
echo ${DTSON[$(($j-1))]}   $ bash array begins at 0
./hindcast_NEST.bash $fatherin $sonin $nestson $j ${DTSON[$(($j-1))]}
#read -p "$*"  #pause        
echo RUNNING CASE $sonin
#mpirun --hostfile /home/fernandotcbarreto/atlasul_operational_noAS/hostfile -np 8 ./romsM $sonin
#./romsS < $sonin
mpirun -np 4 ./romsM_tide $sonin
#read -p "$*"  #pause
done

#./hindcast_NEST.bash $fatherin $sonin $nestson $j
