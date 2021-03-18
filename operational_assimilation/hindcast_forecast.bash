source parameters_operational.bash

BEGINRUN=TRUE

numdaysrun=($(
python - <<EOF
print($numdays - 0.5)
EOF
))

####HINDCAST

nameini=operational_in.in                              
rstday=0.5


inig=-$numdays
endg=0

inim=$inig         # 
endm=$(($endg))  #simulation will start at 12 PM, mercator reference



#################SETUP GFS

#clean GFS data directory
cd $gfsdata
#rm -rf *
cd $maindir
################SETUP MERCATOR

#clean Mercator data directory
cd $mercatordata
#rm -rf *
cd $maindir


# Download the last days of gfs


sed  "s/ini=.*/ini='${inig} days'/" $getgfs | \
sed  "s/tmax=.*/tmax=`date --date "${endg} days" +0z%d%b%Y`/"  > download_gfs.sh

chmod 777 download_gfs.sh
./download_gfs.sh

mv GFS_F* $gfsdata

cd $gfsdata
mv ../grads* .

frctime=`date --date "$inig days" +%Y%m%d`

atmvars=(DLWRAD DSWRAD PAIR QAIR RAIN TAIR USWRAD UWIND VWIND)

for vari in ${atmvars[@]}; do
 filename=GFS_F_${vari}_${frctime}.nc
 while [ !  -f $filename ] || [ ! -s $filename ];do
  echo "missing file $filename, download again"
  grads -lbc "run grads_ncep_$vari.gs"
 done
done                                                           

rm -rf grads*
cd ..

 
##################### Download the last days of Mercator

sed  "s/xmin=.*/xmin=$xmin/" $getmercator | \
sed  "s/xmax=.*/xmax=$xmax/" | \
sed  "s/ymin=.*/ymin=$ymin/" | \
sed  "s/ymax=.*/ymax=$ymax/" | \
sed  "s/for ((.*/for ((  i = ${inim};  i <=${endm} ;  i++  ))/" > download_mercator.sh


chmod 777 download_mercator.sh

./download_mercator.sh

mv MYOCEAN* $mercatordata



################ MAKE ATM BOUNDARY CONDITIONS
cd ${gfsdata}

cp ${foratm}* .

for filename in `ls GFS_F*`; do
 vari=`echo $filename | cut -c7- | rev | cut -c13- | rev`
 echo $vari
 sed -i "0,/${vari}.*/{s/${vari}.*/${vari}file=\'$filename\'/}" ncep2roms_version_1_gfs.py
done

sed -i "0,/input_path=.*/{s@input_path=.*@input_path=\'`pwd`\/\'@}" parameters.py

python ncep2roms_version_1_gfs.py

rm -rf GFS_F*

sed "s@FRCPATH@${gfsdata}@g" ${maindir}/${nameini} |
sed "s/FRCDATE/$frctime/g"> $newini
mv $newini $maindir
cd ..



#############MAKE LATERAL BOUNDARY CONDITIONS

cd $mercatordata

cp ${forocean}* .

cp ../$fathergrid . 

inimdate=`date --date "${inim} days" +%Y-%m-%d`

endmdate=`date --date "${endm} days" +%Y-%m-%d`

boundtime=`date --date "$inim days" +%Y%m%d`


 sed -i "0,/timeini.*/{s/timeini.*/timeini = \'$inimdate 00:00:00\'/}" parameters_bry_in.py

 sed -i "0,/timeend.*/{s/timeend.*/timeend = \'$endmdate 00:00:00\'/}" parameters_bry_in.py
 
 sed -i "0,/input_path=.*/{s@input_path=.*@input_path=\'`pwd`\/MYOCEAN_AZUL_FORECAST_\'@}" parameters_bry_in.py

 sed -i "0,/run_name.*/{s@run_name.*@run_name =\'$boundtime\'@}" parameters_bry_in.py

 sed -i "0,/fname_grd.*/{s@fname_grd.*@fname_grd =\'$fathergrid\'@}" parameters_bry_in.py

 sed -i "0,/rotated=.*/{s@rotated=.*@rotated=False@}" parameters_bry_in.py

 sed -i "0,/cutted=.*/{s@cutted=.*@cutted=True@}" parameters_bry_in.py


python myocean_2_roms_bry_cut.py


TIME_REF=`date --date "$inig days" +%Y%m%d`.5D0   #begin at 12 PM

INI_REFST=`date --date "$inig days" +%Y-%m-%d`
ff=`echo $INI_REFST 12:00:00`

nsegs=($(
python - <<EOF
from matplotlib import dates
print((dates.datestr2num('$ff') - dates.datestr2num('1900-01-01'))*24*60*60)
EOF
))


sed -i "0,/t0 =.*/{s@t0 =.*@t0 =$nsegs@}" myocean_2_roms_ini.py
sed -i "0,/'units', 's')/{s@'units', 's')@'units', 'seconds since 1900-01-01 00:00:00')@}" myocean_2_roms_ini.py

python myocean_2_roms_ini.py

if [ $NUDGECLIM == TRUE ];then
python roms_hycom_intepolation_2_no_stationary_myocean_clm.py
fi

#rm -rf MYOCEAN_AZUL*

sed "s@BRYPATH@${mercatordata}${boundtime}_bry.nc@g" ${maindir}/${newini} |
sed "s@INIPATH@${mercatordata}${boundtime}_ini.nc@g"> ${newini}

if [ $NUDGECLIM == TRUE ];then
sed -i "s@CLMPATH@${mercatordata}${boundtime}_clm.nc@g" ${newini}
sed -i "s@NUDPATH@${dadnud}@g" ${newini}

fi

mv $newini $maindir


cd ..

################################# PROCESSING .IN FILE


#outputing grep to an bash array (using ()

# g=($(grep FRCNAME $INIfile))


#hibrid python and bash

sizegrid=($(
python - <<EOF
from netCDF4 import Dataset
file=Dataset('${fathergrid}')
ntimes=$numdaysrun*24*60*60/($DT)
rsttime=$rstday*24*60*60/($DT)
histime=$hisinterson*60*60/($DT)
#print(some_text)
print(file['lat_rho'][:].shape[0]-2)
print(file['lat_rho'][:].shape[1]-2)
print(ntimes)
print(rsttime)
print(histime)
file
EOF
))

sed -i "0,/GRDNAME ==.*/{s/GRDNAME ==.*/GRDNAME == $fathergrid/}" ${newini}

sed -i "0,/Lm ==.*/{s/Lm ==.*/Lm == ${sizegrid[1]}/}" ${newini}

sed -i "0,/Mm ==.*/{s/Mm ==.*/Mm == ${sizegrid[0]}/}" ${newini}

sed -i "0,/NTIMES ==.*/{s/NTIMES ==.*/NTIMES == ${sizegrid[2]}/}" ${newini}

sed -i "0,/NRST ==.*/{s/NRST ==.*/NRST == ${sizegrid[3]}/}" ${newini}

sed -i "0,/NAVG ==.*/{s/NAVG ==.*/NAVG == 0/}" ${newini}   #no outputing avg

sed -i "0,/DT ==.*/{s/DT ==.*/DT == ${DT}/}" ${newini}


TIME_REF=`date --date "$inig days" +%Y%m%d`.5D0   #begin at 12 PM

INI_REFST=`date --date "$inig days" +%Y-%m-%d`
ff=`echo $INI_REFST 12:00:00`

startd=($(
python - <<EOF
from matplotlib import dates
print(dates.datestr2num('$ff') - dates.datestr2num('1900-01-01'))
EOF
))

TIME_REFSTR=19000101.0D0   #begin at 12 PM

sed -i "0,/DSTART =.*/{s/DSTART =.*/DSTART = $startd/}" ${newini}

sed -i "0,/TIME_REF =.*/{s/TIME_REF =.*/TIME_REF = $TIME_REFSTR/}" ${newini}

sed -i "0,/NHIS ==.*/{s/NHIS ==.*/NHIS == ${sizegrid[4]}/}" ${newini}   

sed -i "0,/NDEFHIS ==.*/{s/NDEFHIS ==.*/NDEFHIS == 0 /}" ${newini}     #only 1 file and no 0001-10 add to his name



sed -i "0,/LcycleRST ==.*/{s/LcycleRST ==.*/LcycleRST == F/}" ${newini}   #need all restarts


##RST AND HIS NAME

TIME_END=`date --date "$endg days" +%Y%m%d`.0D0   #begin at 00 AM

sed -i "0,/HISNAME ==.*/{s/HISNAME ==.*/HISNAME == HIS_FILE_`echo ${TIME_REF} | sed "s/\./_/"`-`echo ${TIME_END} | sed "s/\./_/"`_hind.nc/}" ${newini}   #no outputing avg

sed -i "0,/RSTNAME ==.*/{s/RSTNAME ==.*/RSTNAME == RST_FILE_`echo ${TIME_END} | sed "s/\./_/"`.nc/}" ${newini}   #no outputing avg


if [ $NUDGECLIM == TRUE ];then
sed -i "0,/Lm2CLM ==.*/{s/Lm2CLM ==.*/Lm2CLM == T/}" ${newini}   #only need last restart
sed -i "0,/Lm3CLM ==.*/{s/Lm3CLM ==.*/Lm3CLM == T/}" ${newini}   #only need last restart
sed -i "0,/LtracerCLM ==.*/{s/LtracerCLM ==.*/LtracerCLM == T T/}" ${newini}   #only need last restart
sed -i "0,/LnudgeM2CLM ==.*/{s/LnudgeM2CLM ==.*/LnudgeM2CLM == T/}" ${newini}   #only need last restart
sed -i "0,/LnudgeM3CLM ==.*/{s/LnudgeM3CLM ==.*/LnudgeM3CLM == T/}" ${newini}   #only need last restart
sed -i "0,/LnudgeTCLM ==.*/{s/LnudgeTCLM ==.*/LnudgeTCLM == T T/}" ${newini}   #only need last restart
fi



./romsS	< $newini

if [ $DoNest == TRUE ];then
  ./nesting_run_hind.bash $newini
fi

RST_FROM_HINDCAST=TRUE


################################################### ANALYSIS


source parameters_operational.bash

rstday=1.


inig=$((-$numdaysas -$dly))
endg=$((0 -$dly))

inim=$(($inig -1))        # always 1 minus gfs
endm=$(($endg ))  #simulation will start at 0 AM, mercator reference

################################# PROCESSING .IN FILE

#hibrid python and bash

sizegrid=($(
python - <<EOF
from netCDF4 import Dataset
file=Dataset('${fathergrid}')
ntimes=$numdaysas*24*60*60/($DT)
rsttime=$rstday*24*60*60/($DT)
histime=$hisinterson*60*60/($DT)
#print(some_text)
print(file['lat_rho'][:].shape[0]-2)
print(file['lat_rho'][:].shape[1]-2)
print(ntimes)
print(rsttime)
print(histime)
file
EOF
))

sed -i "0,/GRDNAME ==.*/{s/GRDNAME ==.*/GRDNAME == $fathergrid/}" ${newini}

sed -i "0,/Lm ==.*/{s/Lm ==.*/Lm == ${sizegrid[1]}/}" ${newini}

sed -i "0,/Mm ==.*/{s/Mm ==.*/Mm == ${sizegrid[0]}/}" ${newini}

sed -i "0,/NTIMES ==.*/{s/NTIMES ==.*/NTIMES == ${sizegrid[2]}/}" ${newini}

sed -i "0,/NRST ==.*/{s/NRST ==.*/NRST == ${sizegrid[3]}/}" ${newini}

sed -i "0,/NAVG ==.*/{s/NAVG ==.*/NAVG == 2000000000/}" ${newini}   #no outputing avg

sed -i "0,/DT ==.*/{s/DT ==.*/DT == ${DT}/}" ${newini}


TIME_REF=`date --date "$inig days" +%Y%m%d`.0D0   #begin at 00 AM

INI_REFST=`date --date "$inig days" +%Y-%m-%d`
ff=`echo $INI_REFST 00:00:00`

startd=($(
python - <<EOF
from matplotlib import dates
print(dates.datestr2num('$ff') - dates.datestr2num('1900-01-01'))
EOF
))

TIME_REFSTR=19000101.0D0   #begin at 12 PM

sed -i "0,/DSTART =.*/{s/DSTART =.*/DSTART = $startd/}" ${newini}

sed -i "0,/TIME_REF =.*/{s/TIME_REF =.*/TIME_REF = $TIME_REFSTR/}" ${newini}

sed -i "0,/NHIS ==.*/{s/NHIS ==.*/NHIS == ${sizegrid[4]}/}" ${newini}   

sed -i "0,/NDEFHIS ==.*/{s/NDEFHIS ==.*/NDEFHIS == 0 /}" ${newini}     #only 1 file and no 0001-10 add to his name


sed -i "0,/LcycleRST ==.*/{s/LcycleRST ==.*/LcycleRST == F/}" ${newini}      #Need all restarts



##RST AND HIS NAME

TIME_END=`date --date "$endg days" +%Y%m%d`.0D0   #begin at 00 AM

inim=$(($endg +dly)) # RST of tomorrow is today plus 1, write new RST NAME

TIME_FOR=`date --date "$inim days" +%Y%m%d`.0D0   #begin at 00 AM

sed -i "0,/HISNAME ==.*/{s/HISNAME ==.*/HISNAME == HIS_FILE_`echo ${TIME_REF} | sed "s/\./_/"`-`echo ${TIME_END} | sed "s/\./_/"`_fore.nc/}" ${newini}   #no outputing avg

sed -i "0,/RSTNAME ==.*/{s/RSTNAME ==.*/RSTNAME == RST_FILE_`echo ${TIME_END} | sed "s/\./_/"`_AS.nc/}" ${newini}   # RST of forecast is today plus one



#####INI FILE COMES FROM RESTART

sed -i "0,/ININAME.*/{s@ININAME.*@ININAME == RST_FILE_`echo ${TIME_FOR} | sed "s/\./_/"`.nc@}"  ${newini}

rstfile=\'RST_FILE_`echo ${TIME_FOR} | sed "s/\./_/"`.nc\'



rstindex=($(
python - <<EOF
from netCDF4 import Dataset
from matplotlib import dates
import numpy as np
file=Dataset($rstfile)
#print(some_text)
print(int((np.where(file['ocean_time'][:]/(24*60*60)==dates.datestr2num('$ff') - dates.datestr2num('1900-01-01'))[0])))
file
EOF
))

rstindex=$((${rstindex[0]}+1))


sed -i "0,/NRREC ==.*/{s/NRREC ==.*/NRREC == $rstindex/}" ${newini}   #no outputing avg


if [ $NUDGECLIM == TRUE ];then
sed -i "0,/Lm2CLM ==.*/{s/Lm2CLM ==.*/Lm2CLM == T/}" ${newini}   #only need last restart
sed -i "0,/Lm3CLM ==.*/{s/Lm3CLM ==.*/Lm3CLM == T/}" ${newini}   #only need last restart
sed -i "0,/LtracerCLM ==.*/{s/LtracerCLM ==.*/LtracerCLM == T T/}" ${newini}   #only need last restart
sed -i "0,/LnudgeM2CLM ==.*/{s/LnudgeM2CLM ==.*/LnudgeM2CLM == T/}" ${newini}   #only need last restart
sed -i "0,/LnudgeM3CLM ==.*/{s/LnudgeM3CLM ==.*/LnudgeM3CLM == T/}" ${newini}   #only need last restart
sed -i "0,/LnudgeTCLM ==.*/{s/LnudgeTCLM ==.*/LnudgeTCLM == T T/}" ${newini}   #only need last restart
fi


./romsS	< $newini

rm *outer0.nc
rnm=$(grep HISNAME $newini | head -n 1 | tr -s ' '| sed "s/' '/''/g" | cut -d '=' -f3 | cut -d '.' -f1)
mv BRSE_fwd_outer1.nc ${rnm}.nc

if [ $DoNest == TRUE ];then
  ./nesting_run_fore_ass.bash $newini
fi


RST_FROM_HINDCAST=FALSE



##############################################################################################
########################## FORECAST


source parameters_operational.bash

nameini=operational_in.in                              
rstday=0.5


inig=$((0 -$dly))
endg=$(($numdaysfor -$dly))

inim=$(($inig -1))         # always 1 minus gfs
endm=$(($endg))  #simulation will start at 12 PM, mercator reference



#################SETUP GFS

#clean GFS data directory
cd $gfsdata
#rm -rf *
cd $maindir
################SETUP MERCATOR

#clean Mercator data directory
cd $mercatordata
#rm -rf *
cd $maindir

# Download the last days of gfs


sed  "s/ini=.*/ini='${inig} days'/" $getgfs | \
sed  "s/tmax=.*/tmax=`date --date "${endg} days" +0z%d%b%Y`/"  > download_gfs.sh

chmod 777 download_gfs.sh
./download_gfs.sh

mv GFS_F* $gfsdata

cd $gfsdata
mv ../grads* .

frctime=`date --date "$inig days" +%Y%m%d`

atmvars=(DLWRAD DSWRAD PAIR QAIR RAIN TAIR USWRAD UWIND VWIND)

for vari in ${atmvars[@]}; do
 filename=GFS_F_${vari}_${frctime}.nc
 while [ !  -f $filename ] || [ ! -s $filename ];do
  echo "missing file $filename, download again"
  grads -lbc "run grads_ncep_$vari.gs"
 done
done                                                           

rm -rf grads*
cd ..

 
##################### Download the last days of Mercator

sed  "s/xmin=.*/xmin=$xmin/" $getmercator | \
sed  "s/xmax=.*/xmax=$xmax/" | \
sed  "s/ymin=.*/ymin=$ymin/" | \
sed  "s/ymax=.*/ymax=$ymax/" | \
sed  "s/for ((.*/for ((  i = ${inim};  i <=${endm} ;  i++  ))/" > download_mercator.sh


chmod 777 download_mercator.sh

./download_mercator.sh

mv MYOCEAN* $mercatordata



################ MAKE ATM BOUNDARY CONDITIONS
cd ${gfsdata}

cp ${foratm}* .

for filename in `ls GFS_F*`; do
 vari=`echo $filename | cut -c7- | rev | cut -c13- | rev`
 echo $vari
 sed -i "0,/${vari}.*/{s/${vari}.*/${vari}file=\'$filename\'/}" ncep2roms_version_1_gfs.py
done

sed -i "0,/input_path=.*/{s@input_path=.*@input_path=\'`pwd`\/\'@}" parameters.py

python ncep2roms_version_1_gfs.py

rm -rf GFS_F*


sed "s@FRCPATH@${gfsdata}@g" ${maindir}/${nameini} |
sed "s/FRCDATE/$frctime/g"> $newini
mv $newini $maindir
cd ..



#############MAKE LATERAL BOUNDARY CONDITIONS

cd $mercatordata

cp ${forocean}* .


inimdate=`date --date "${inim} days" +%Y-%m-%d`

endmdate=`date --date "${endm} days" +%Y-%m-%d`

boundtime=`date --date "$inim days" +%Y%m%d`


 sed -i "0,/timeini.*/{s/timeini.*/timeini = \'$inimdate 00:00:00\'/}" parameters_bry_in.py

 sed -i "0,/timeend.*/{s/timeend.*/timeend = \'$endmdate 00:00:00\'/}" parameters_bry_in.py
 
 sed -i "0,/input_path=.*/{s@input_path=.*@input_path=\'`pwd`\/MYOCEAN_AZUL_FORECAST_\'@}" parameters_bry_in.py

 sed -i "0,/run_name.*/{s@run_name.*@run_name =\'$boundtime\'@}" parameters_bry_in.py

 sed -i "0,/fname_grd.*/{s@fname_grd.*@fname_grd =\'$fathergrid\'@}" parameters_bry_in.py

 sed -i "0,/rotated=.*/{s@rotated=.*@rotated=False@}" parameters_bry_in.py

 sed -i "0,/cutted=.*/{s@cutted=.*@cutted=True@}" parameters_bry_in.py


python myocean_2_roms_bry_cut.py

if [ $NUDGECLIM == TRUE ];then
python roms_hycom_intepolation_2_no_stationary_myocean_clm.py
fi

#rm -rf MYOCEAN_AZUL*


sed "s@BRYPATH@${mercatordata}${boundtime}_bry.nc@g" ${maindir}/${newini} > ${newini}

if [ $NUDGECLIM == TRUE ];then
sed -i "s@CLMPATH@${mercatordata}${boundtime}_clm.nc@g" ${newini}
sed -i "s@NUDPATH@${dadnud}@g" ${newini}

fi


mv $newini $maindir


cd ..

################################# PROCESSING .IN FILE


#outputing grep to an bash array (using ()

# g=($(grep FRCNAME $INIfile))


#hibrid python and bash

sizegrid=($(
python - <<EOF
from netCDF4 import Dataset
file=Dataset('${fathergrid}')
ntimes=$numdaysfor*24*60*60/($DT)
rsttime=$rstday*24*60*60/($DT)
histime=$hisinterson*60*60/($DT)
#print(some_text)
print(file['lat_rho'][:].shape[0]-2)
print(file['lat_rho'][:].shape[1]-2)
print(ntimes)
print(rsttime)
print(histime)
file
EOF
))

sed -i "0,/GRDNAME ==.*/{s/GRDNAME ==.*/GRDNAME == $fathergrid/}" ${newini}

sed -i "0,/Lm ==.*/{s/Lm ==.*/Lm == ${sizegrid[1]}/}" ${newini}

sed -i "0,/Mm ==.*/{s/Mm ==.*/Mm == ${sizegrid[0]}/}" ${newini}

sed -i "0,/NTIMES ==.*/{s/NTIMES ==.*/NTIMES == ${sizegrid[2]}/}" ${newini}

#sed -i "0,/NRST ==.*/{s/NRST ==.*/NRST == ${sizegrid[3]}/}" ${newini}
sed -i "0,/NRST ==.*/{s/NRST ==.*/NRST == 0/}" ${newini}

sed -i "0,/NAVG ==.*/{s/NAVG ==.*/NAVG == 0/}" ${newini}   #no outputing avg

sed -i "0,/DT ==.*/{s/DT ==.*/DT == ${DT}/}" ${newini}


TIME_REF=`date --date "$inig days" +%Y%m%d`.0D0   #begin at 00 AM

INI_REFST=`date --date "$inig days" +%Y-%m-%d`
ff=`echo $INI_REFST 00:00:00`

startd=($(
python - <<EOF
from matplotlib import dates
print(dates.datestr2num('$ff') - dates.datestr2num('1900-01-01'))
EOF
))

TIME_REFSTR=19000101.0D0   #begin at 12 PM

 
sed -i "0,/DSTART =.*/{s/DSTART =.*/DSTART = $startd/}" ${newini}

sed -i "0,/TIME_REF =.*/{s/TIME_REF =.*/TIME_REF = $TIME_REFSTR/}" ${newini}

sed -i "0,/NHIS ==.*/{s/NHIS ==.*/NHIS == ${sizegrid[4]}/}" ${newini}   

sed -i "0,/NDEFHIS ==.*/{s/NDEFHIS ==.*/NDEFHIS == 0 /}" ${newini}     #only 1 file and no 0001-10 add to his name


sed -i "0,/LcycleRST ==.*/{s/LcycleRST ==.*/LcycleRST == F/}" ${newini}      #Need all restarts



#####INI FILE COMES FROM RESTART

sed -i "s@INIPATH@${maindir}/RST_FILE_`echo ${TIME_REF} | sed "s/\./_/"`_AS.nc@g" ${newini}
rstfile=\'RST_FILE_`echo ${TIME_REF} | sed "s/\./_/"`_AS.nc\'



##RST AND HIS NAME

TIME_END=`date --date "$endg days" +%Y%m%d`.0D0   #begin at 00 AM

inim=$(($inig + 1)) # RST of tomorrow is today plus 1, write new RST NAME

TIME_FOR=`date --date "$inim days" +%Y%m%d`.0D0   #begin at 00 AM

sed -i "0,/HISNAME ==.*/{s/HISNAME ==.*/HISNAME == HIS_FILE_`echo ${TIME_REF} | sed "s/\./_/"`-`echo ${TIME_END} | sed "s/\./_/"`_fore.nc/}" ${newini}   #no outputing avg

sed -i "0,/RSTNAME ==.*/{s/RSTNAME ==.*/RSTNAME == RST_FILE_`echo ${TIME_FOR} | sed "s/\./_/"`.nc/}" ${newini}   # RST of forecast is today plus one



rstindex=($(
python - <<EOF
from netCDF4 import Dataset
import numpy as np
from matplotlib import dates
file=Dataset($rstfile)
#print(some_text)
print(int((np.where(file['ocean_time'][:]/(24*60*60)==dates.datestr2num('$ff') - dates.datestr2num('1900-01-01'))[0])))
file
EOF
))


rstindex=$((${rstindex[0]}+1))


sed -i "0,/NRREC ==.*/{s/NRREC ==.*/NRREC == $rstindex/}" ${newini}   #no outputing avg


if [ $NUDGECLIM == TRUE ];then
sed -i "0,/Lm2CLM ==.*/{s/Lm2CLM ==.*/Lm2CLM == T/}" ${newini}   #only need last restart
sed -i "0,/Lm3CLM ==.*/{s/Lm3CLM ==.*/Lm3CLM == T/}" ${newini}   #only need last restart
sed -i "0,/LtracerCLM ==.*/{s/LtracerCLM ==.*/LtracerCLM == T T/}" ${newini}   #only need last restart
sed -i "0,/LnudgeM2CLM ==.*/{s/LnudgeM2CLM ==.*/LnudgeM2CLM == T/}" ${newini}   #only need last restart
sed -i "0,/LnudgeM3CLM ==.*/{s/LnudgeM3CLM ==.*/LnudgeM3CLM == T/}" ${newini}   #only need last restart
sed -i "0,/LnudgeTCLM ==.*/{s/LnudgeTCLM ==.*/LnudgeTCLM == T T/}" ${newini}   #only need last restart
fi


./romsS	< $newini

if [ $DoNest == TRUE ];then
  ./nesting_run_fore.bash $newini
fi


RST_FROM_HINDCAST=FALSE








#a=`date --date "0 days" +%H`
#if  [ $a -gt 10 ] && [ $a -lt 11 ] 
#then 
#echo 1432;
#fi

#NTIMES

#echo ${sizegrid[@]}.0D0


#sed -i "0,/MY_PROJECT_DIR.*/{s@MY_PROJECT_DIR.*@MY_PROJECT_DIR `pwd`@}" build_roms.sh



# grep RSTNAME ocean_in_forecast.in | head -n 1 |  tr -s ' '| cut -d ' ' -f4   ## tr -s ' ' get rid of blank spaces
 
# d=2015-01-01
# while [ "$d" != 2015-02-20 ]; do 
#   echo $d
#   d=$(date -I -d "$d + 1 day")
#done