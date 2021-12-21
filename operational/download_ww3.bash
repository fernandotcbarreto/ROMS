#!/bin/bash -x

xmin=-55
xmax=-30
ymin=-40
ymax=-15

username='fbarreto'
password='Mg05031990*'

        ano=`date --date 'today' +%Y`
        mes=`date --date 'today' +%m`
        dia=`date --date 'today' +%d`

        anomin=`date --date "+1 days" +%Y`
        mesmin=`date --date "+1 days" +%m`
        diamin=`date --date "+1 days" +%d`
        anomax=`date --date "+2 days" +%Y`
        mesmax=`date --date "+2 days" +%m`
        diamax=`date --date "+2 days" +%d`
        
        
        filename=previsao_preliminar_$anomin$mesmin$diamin.nc
        
        while [ !  -f $filename ] || [ ! -s $filename ];do

          echo DOWNLOADING $filename
          
	  #python -m motuclient --motu http://nrt.cmems-du.eu/motu-web/Motu --service-id GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS --product-id global-analysis-forecast-phy-001-024 --longitude-min $xmin --longitude-max $xmax --latitude-min $ymin --latitude-max $ymax --date-min "$anomin-$mesmin-$diamin 12:00:00" --date-max "$anomax-$mesmax-$diamax 12:00:00" --depth-min 0.493 --depth-max 5727.9169921875 --variable thetao --variable so --variable zos --variable uo --variable vo --out-dir `pwd` --out-name MYOCEAN_AZUL_FORECAST_$anomin$mesmin$diamin.nc --user $username --pwd $password

      python -m motuclient --motu https://nrt.cmems-du.eu/motu-web/Motu --service-id GLOBAL_ANALYSIS_FORECAST_WAV_001_027-TDS --product-id global-analysis-forecast-wav-001-027 --longitude-min $xmin --longitude-max $xmax --latitude-min $ymin --latitude-max $ymax --date-min "$anomin-$mesmin-$diamin 00:00:00" --date-max "$anomax-$mesmax-$diamax 00:00:00" --variable VHM0 --variable VSDX --variable VSDY --out-dir `pwd` --out-name previsao_preliminar_$anomin$mesmin$diamin.nc --user $username --pwd $password
       done



