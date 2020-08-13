#!/bin/sh -x

# caminho da sua pasta
#cd /data/scripts/MORROSQUILLO/
# define dia inicial e final da previsao

tmin=`date --date 'Today' +0z%d%b%Y`
tmax=`date --date '+9 days' +0z%d%b%Y`
tstring=`date --date 'Today' +%Y%m%d`

####10 gives an error

# Definindo latitude (indice)

iymin=230 #
iymax=310 #

# Definindo longitude (indice)

ixmin=1220 #
ixmax=1350 #


# Cria arquivo gs da componente U para executar no GrADSDAP
#
cat > grads_ncep_UWIND.gs << EOF
* Abre o opendap

'sdfopen https://nomads.ncep.noaa.gov/dods/gfs_0p25/gfs$tstring/gfs_0p25_00z'    


* Define o nome do arquivo binario que sera gerado pelo GrADS

'set sdfwrite GFS_F_UWIND_$tstring.nc'


*'set y 'iymin' 'iymax
'set y $iymin $iymax'

* 'set x 'ixmin' 'ixmax   

'set x $ixmin $ixmax'

* Definindo tempo

'set time $tmin $tmax'

'define ugrd10m = ugrd10m'

'sdfwrite ugrd10m'


*Fecha arquivos
'quit'

EOF
#

# Cria arquivo gs da componente V para executar no GrADSDAP

cat > grads_ncep_VWIND.gs << EOF
* Abre o opendap
    
'sdfopen https://nomads.ncep.noaa.gov/dods/gfs_0p25/gfs$tstring/gfs_0p25_00z'   


* Define o nome do arquivo binario que sera gerado pelo GrADS

'set sdfwrite  GFS_F_VWIND_$tstring.nc'


*'set y 'iymin' 'iymax
'set y $iymin $iymax'

* 'set x 'ixmin' 'ixmax   

'set x $ixmin $ixmax'

* Definindo tempo

'set time $tmin $tmax'

'define vgrd10m = vgrd10m'

'sdfwrite vgrd10m'


*Fecha arquivos
'quit'

EOF

#

# Cria arquivo gs de precipitacao executar no GrADSDAP

cat > grads_ncep_RAIN.gs << EOF
* Abre o opendap
    
'sdfopen https://nomads.ncep.noaa.gov/dods/gfs_0p25/gfs$tstring/gfs_0p25_00z'   

* Define o nome do arquivo binario que sera gerado pelo GrADS

'set sdfwrite  GFS_F_RAIN_$tstring.nc'


*'set y 'iymin' 'iymax
'set y $iymin $iymax'

* 'set x 'ixmin' 'ixmax   

'set x $ixmin $ixmax'

* Definindo tempo

'set time $tmin $tmax'

'define pratesfc = pratesfc'

'sdfwrite pratesfc'


*Fecha arquivos
'quit'

EOF

#

# Cria arquivo gs de down swave para executar no GrADSDAP

cat > grads_ncep_DSWRAD.gs << EOF
* Abre o opendap
    
'sdfopen https://nomads.ncep.noaa.gov/dods/gfs_0p25/gfs$tstring/gfs_0p25_00z'   


* Define o nome do arquivo binario que sera gerado pelo GrADS

'set sdfwrite  GFS_F_DSWRAD_$tstring.nc'


*'set y 'iymin' 'iymax
'set y $iymin $iymax'

* 'set x 'ixmin' 'ixmax   

'set x $ixmin $ixmax'

* Definindo tempo

'set time $tmin $tmax'

'define dswrfsfc = dswrfsfc'

'sdfwrite dswrfsfc'


*Fecha arquivos
'quit'

EOF

#

# Cria arquivo gs de down lwave para executar no GrADSDAP

cat > grads_ncep_DLWRAD.gs << EOF
* Abre o opendap
    
'sdfopen https://nomads.ncep.noaa.gov/dods/gfs_0p25/gfs$tstring/gfs_0p25_00z'   

* Define o nome do arquivo binario que sera gerado pelo GrADS

'set sdfwrite GFS_F_DLWRAD_$tstring.nc'


*'set y 'iymin' 'iymax
'set y $iymin $iymax'

* 'set x 'ixmin' 'ixmax   

'set x $ixmin $ixmax'

* Definindo tempo

'set time $tmin $tmax'

'define dlwrfsfc = dlwrfsfc'

'sdfwrite dlwrfsfc'


*Fecha arquivos
'quit'

EOF

#

# Cria arquivo gs de pressao para executar no GrADSDAP

cat > grads_ncep_PAIR.gs << EOF
* Abre o opendap
    
'sdfopen https://nomads.ncep.noaa.gov/dods/gfs_0p25/gfs$tstring/gfs_0p25_00z'   


* Define o nome do arquivo binario que sera gerado pelo GrADS

'set sdfwrite GFS_F_PAIR_$tstring.nc'


*'set y 'iymin' 'iymax
'set y $iymin $iymax'

* 'set x 'ixmin' 'ixmax   

'set x $ixmin $ixmax'

* Definindo tempo

'set time $tmin $tmax'

'define pressfc = pressfc'

'sdfwrite pressfc'


*Fecha arquivos
'quit'

EOF

#

# Cria arquivo gs de umidade para executar no GrADSDAP

cat > grads_ncep_QAIR.gs << EOF
* Abre o opendap
    
'sdfopen https://nomads.ncep.noaa.gov/dods/gfs_0p25/gfs$tstring/gfs_0p25_00z'   


* Define o nome do arquivo binario que sera gerado pelo GrADS

'set sdfwrite GFS_F_QAIR_$tstring.nc'


*'set y 'iymin' 'iymax
'set y $iymin $iymax'

* 'set x 'ixmin' 'ixmax   

'set x $ixmin $ixmax'

* Definindo tempo

'set time $tmin $tmax'

'define rh2m = rh2m'

'sdfwrite rh2m'


*Fecha arquivos
'quit'

EOF

#

# Cria arquivo gs de temperatura para executar no GrADSDAP

cat > grads_ncep_TAIR.gs << EOF
* Abre o opendap
    
'sdfopen https://nomads.ncep.noaa.gov/dods/gfs_0p25/gfs$tstring/gfs_0p25_00z'   


* Define o nome do arquivo binario que sera gerado pelo GrADS

'set sdfwrite GFS_F_TAIR_$tstring.nc'


*'set y 'iymin' 'iymax
'set y $iymin $iymax'

* 'set x 'ixmin' 'ixmax   

'set x $ixmin $ixmax'

* Definindo tempo

'set time $tmin $tmax'

'define tmp2m = tmp2m'

'sdfwrite tmp2m'


*Fecha arquivos
'quit'

EOF

#

# Cria arquivo gs de up swave para executar no GrADSDAP

cat > grads_ncep_USWRAD.gs << EOF
* Abre o opendap
    
'sdfopen https://nomads.ncep.noaa.gov/dods/gfs_0p25/gfs$tstring/gfs_0p25_00z'   


* Define o nome do arquivo binario que sera gerado pelo GrADS

'set sdfwrite GFS_F_USWRAD_$tstring.nc'


*'set y 'iymin' 'iymax
'set y $iymin $iymax'

* 'set x 'ixmin' 'ixmax   

'set x $ixmin $ixmax'

* Definindo tempo

'set time $tmin $tmax'

'define uswrfsfc = uswrfsfc'

'sdfwrite uswrfsfc'


*Fecha arquivos
'quit'

EOF

grads -lbc "run grads_ncep_UWIND.gs"
grads -lbc "run grads_ncep_VWIND.gs"
grads -lbc "run grads_ncep_RAIN.gs"
grads -lbc "run grads_ncep_DSWRAD.gs"
grads -lbc "run grads_ncep_DLWRAD.gs"
grads -lbc "run grads_ncep_PAIR.gs"
grads -lbc "run grads_ncep_QAIR.gs"
grads -lbc "run grads_ncep_TAIR.gs"
grads -lbc "run grads_ncep_USWRAD.gs"

# VERIFICACAO
echo '--------------------------------'
echo 'INICIANDO VERIFICACAO DE ARQUIVO'
echo '--------------------------------'
perf_size=234900
files=GFS_F_*.nc
for var in $files; do
        size=$(stat -c%s "$var")
        vari=`echo $var | cut -c7- | rev | cut -c14- | rev`
        echo $vari $size
        while [ $size -lt $perf_size ]; do
                echo 'Arquivo '$var' avariado, baixando de novo...'
                grads -lbc "run grads_ncep_$vari.gs"
                size=$(stat -c%s "$var")
        done
done

echo '--------------------------------'
echo ' FIM DA VERIFICACAO DE ARQUIVO  '
echo '--------------------------------'


mv GFS_F* /rackstation2/dados/Modelos/MORROSQUILLO/GFS/

#rm GFS_F*
