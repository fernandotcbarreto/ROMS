source parameters_operational.bash

fatherin=$1
nestin1=$2
nestson1=$3
ind=$4
DTSON=$5

echo $fatherin
echo $nestson1

DT=$(grep DT $fatherin | head -n 1 | tr -s ' '| sed "s/' '/''/g" | cut -d '=' -f3)

echo $DT

NTIMES=$(grep NTIMES $fatherin | head -n 1 | tr -s ' '| sed "s/' '/''/g" | cut -d '=' -f3)

NRST=$(grep NRST $fatherin | head -n 1 | tr -s ' '| sed "s/' '/''/g" | cut -d '=' -f3)



####

#############MAKE LATERAL BOUNDARY CONDITIONS

cp $nestson1 $mercatordata

HISNAME=$(grep HISNAME $fatherin | head -n 1 | tr -s ' '| sed "s/' '/''/g" | cut -d '=' -f3 | cut -d '.' -f1)
RSTNAME=$(grep RSTNAME $fatherin | head -n 1 | tr -s ' '| sed "s/' '/''/g" | cut -d '=' -f3 | cut -d '.' -f1)

TIME_REF=$(grep TIME_REF $fatherin | head -n 1 | tr -s ' ' | cut -d '=' -f2)
echo $TIME_REF

#echo $HISNAME
#echo $RSTNAME

cp ${HISNAME}.nc $mercatordata

cd $mercatordata

cp ${forocean}* .

sed -i "0,/coup_files.*/{s@coup_files.*@coup_files =\[\'`echo ${HISNAME}`.nc\'\]@}" parameters_bry_in.py

sed -i "0,/fname_grd_son.*/{s@fname_grd_son.*@fname_grd_son =\'`echo ${nestson1}`\'@}" parameters_bry_in.py

sed -i "0,/run_name_son.*/{s@run_name_son.*@run_name_son =\'`echo ${TIME_REF} | cut -d '.' -f1`_NEST_${ind}\'@}" parameters_bry_in.py

sed -i "0,/rotated=.*/{s@rotated=.*@rotated=$rotatedid@}" parameters_bry_in.py

sed -i "0,/cutted=.*/{s@cutted=.*@cutted=$cuttedid@}" parameters_bry_in.py

sed -i "0,/lim=.*/{s@lim=.*@lim=$lim@}" parameters_bry_in.py


python roms_2_roms_bry_cut.py
python roms_2_roms_ini.py

sed "0,/BRYNAME.*/{s@BRYNAME.*@BRYNAME == ${mercatordata}`echo ${TIME_REF} | cut -d '.' -f1`_NEST_${ind}_bry.nc@}" ${maindir}/${fatherin} |

sed "0,/ININAME.*/{s@ININAME.*@ININAME == ${mercatordata}`echo ${TIME_REF} | cut -d '.' -f1`_NEST_${ind}_ini.nc@}" > ${nestin1}

if [ $NUDGECLIM == TRUE ];then
  python roms_2_roms_clm.py
  sed -i "0,/CLMNAME.*/{s@CLMNAME.*@CLMNAME == ${mercatordata}`echo ${TIME_REF} | cut -d '.' -f1`_NEST_${ind}_clm.nc@}" ${nestin1}
  var="son_nud_${ind}"
  sed -i "0,/NUDNAME.*/{s@NUDNAME.*@NUDNAME == ${!var} @}" ${nestin1}
fi

var="visc2_${ind}"
sed -i "0,/VISC2.*/{s@VISC2.*@VISC2 == ${!var}.0d0 @}" ${nestin1}

mv $nestin1 $maindir

cd ..


####################################################
echo $NTIMES

echo $NRST

timeval=($(
python - <<EOF
import numpy as np
ndays=(float(${NTIMES})*float(${DT}))/(24.*60.*60.)
ndaysrst=(float(${NRST})*float(${DT}))/(24.*60.*60.)
print(ndays, ndaysrst)
EOF
))

numdays=${timeval[0]}
rstday=${timeval[1]}

echo $numdays
echo $rstday

sizegrid=($(
python - <<EOF
from netCDF4 import Dataset
file=Dataset('$nestson1')
ntimes=$numdays*24*60*60/($DTSON)
rsttime=$rstday*24*60*60/($DTSON)
histime=$hisinterson*60*60/($DTSON)
#print(some_text)
print(file['lat_rho'][:].shape[0]-2)
print(file['lat_rho'][:].shape[1]-2)
print(ntimes)
print(rsttime)
print(histime)
file
EOF
))

echo ${sizegrid[@]}

sed -i "0,/Lm ==.*/{s/Lm ==.*/Lm == ${sizegrid[1]}/}" ${nestin1}

sed -i "0,/Mm ==.*/{s/Mm ==.*/Mm == ${sizegrid[0]}/}" ${nestin1}

sed -i "0,/NTIMES ==.*/{s/NTIMES ==.*/NTIMES == ${sizegrid[2]}/}" ${nestin1}

sed -i "0,/NRST ==.*/{s/NRST ==.*/NRST == ${sizegrid[3]}/}" ${nestin1}

sed -i "0,/NHIS ==.*/{s/NHIS ==.*/NHIS == ${sizegrid[4]}/}" ${nestin1}   #no outputing avg

sed -i "0,/NAVG ==.*/{s/NAVG ==.*/NAVG == 0/}" ${nestin1}   #no outputing avg

sed -i "0,/NDEFHIS ==.*/{s/NDEFHIS ==.*/NDEFHIS == 0/}" ${nestin1}   #no outputing avg

sed -i "0,/NDIA ==.*/{s/NDIA ==.*/NDIA == 0/}" ${nestin1}   #no outputing avg

sed -i "0,/DT ==.*/{s/DT ==.*/DT == ${DTSON}/}" ${nestin1}

sed -i "0,/GRDNAME ==.*/{s/GRDNAME ==.*/GRDNAME == `echo $nestson1`/}" ${nestin1}   #no outputing avg
sed -i "0,/HISNAME ==.*/{s/HISNAME ==.*/HISNAME == `echo $HISNAME`_NEST_${ind}.nc/}" ${nestin1}   #no outputing avg
sed -i "0,/RSTNAME ==.*/{s/RSTNAME ==.*/RSTNAME == `echo $RSTNAME`_NEST_${ind}.nc/}" ${nestin1}   #no outputing avg

if [ $NUDGECLIM == FALSE ];then
  sed -i "0,/Lm2CLM ==.*/{s/Lm2CLM ==.*/Lm2CLM == F/}" ${nestin1}   #only need last restart
  sed -i "0,/Lm3CLM ==.*/{s/Lm3CLM ==.*/Lm3CLM == F/}" ${nestin1}   #only need last restart
  sed -i "0,/LtracerCLM ==.*/{s/LtracerCLM ==.*/LtracerCLM == F F/}" ${nestin1}   #only need last restart
  sed -i "0,/LnudgeM2CLM ==.*/{s/LnudgeM2CLM ==.*/LnudgeM2CLM == F/}" ${nestin1}   #only need last restart
  sed -i "0,/LnudgeM3CLM ==.*/{s/LnudgeM3CLM ==.*/LnudgeM3CLM == F/}" ${nestin1}   #only need last restart
  sed -i "0,/LnudgeTCLM ==.*/{s/LnudgeTCLM ==.*/LnudgeTCLM == F F/}" ${nestin1}   #only need last restart
fi


