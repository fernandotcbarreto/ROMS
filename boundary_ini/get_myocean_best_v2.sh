#!/bin/bash -x

xmin=-55
xmax=-30
ymin=-40
ymax=-15

username='fbarreto'
password='Mg05031990*'
for ((  i = -7;  i <=0 ;  i++  ))
do
        ano=`date --date 'today' +%Y`
        mes=`date --date 'today' +%m`
        dia=`date --date 'today' +%d`

        anomin=`date --date "+$i days" +%Y`
        mesmin=`date --date "+$i days" +%m`
        diamin=`date --date "+$i days" +%d`
        anomax=`date --date "+$i days" +%Y`
        mesmax=`date --date "+$i days" +%m`
        diamax=`date --date "+$i days" +%d`
        
        
        filename=MYOCEAN_AZUL_FORECAST_$anomin$mesmin$diamin.nc
        
        while [ !  -f $filename ] || [ ! -s $filename ];do

          echo DOWNLOADING $filename
          
	  python -m motuclient --motu http://nrt.cmems-du.eu/motu-web/Motu --service-id GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS --product-id global-analysis-forecast-phy-001-024 --longitude-min $xmin --longitude-max $xmax --latitude-min $ymin --latitude-max $ymax --date-min "$anomin-$mesmin-$diamin 12:00:00" --date-max "$anomax-$mesmax-$diamax 12:00:00" --depth-min 0.493 --depth-max 5727.9169921875 --variable thetao --variable so --variable zos --variable uo --variable vo --out-dir `pwd` --out-name MYOCEAN_AZUL_FORECAST_$anomin$mesmin$diamin.nc --user $username --pwd $password
       
        done


done


# EU UTILIZEI 

#python -m motuclient --motu http://nrt.cmems-du.eu/motu-web/Motu --service-id GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS --product-id global-analysis-forecast-phy-001-024 --longitude-min -50 --longitude-max -35 --latitude-min -30 --latitude-max -10 --date-min "2020-04-15 12:00:00" --date-max "2020-04-15 12:00:00" --depth-min 0.493 --depth-max 5727.918000000001 --variable thetao --variable so --variable zos --variable uo --variable vo --out-dir /home/fernando --out-name saida_january.nc --user fbarreto --pwd Mg05031990*


#python -m motuclient --motu http://nrt.cmems-du.eu/motu-web/Motu --service-id GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS --product-id global-analysis-forecast-phy-001-024 --longitude-min -50 --longitude-max -35 --latitude-min -30 --latitude-max -10 --date-min "2020-01-15 12:00:00" --date-max "2020-01-15 12:00:00" --depth-min 0.493 --depth-max 5727.918000000001 --variable thetao --variable so --variable zos --variable uo --variable vo --out-dir /home/fernando --out-name saida_january.nc --user fbarreto --pwd Mg05031990*


