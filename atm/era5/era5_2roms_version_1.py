from parameters import *
from datetime import date
from netCDF4 import Dataset as dat
import numpy as np
from forclass import ncgene
import matplotlib.dates as dates
import matplotlib.pyplot as plt

UWINDfile='download_u_10_2016.nc'
VWINDfile='download_v_10_2016.nc'
RAINfile='precipitation_2016.nc'
DLWRADfile='downward_long_2016.nc'
TAIRfile='download_temperature_2016.nc'
PAIRfile='download_p_2016.nc'
DSWRADfile='net_solar_2016.nc'
QAIRfile='download_humidity_2016.nc'

vari='2017'

fnin=dat(input_path + UWINDfile)                  #name of U wind file
uwindf=np.squeeze(fnin["u10"][::])
uwindf=np.transpose(uwindf,(2,1,0))


fninv=dat(input_path + VWINDfile)                  #name of V wind file
vwindf=np.squeeze(fninv["v10"][::])
vwindf=np.transpose(vwindf,(2,1,0))


fninr=dat(input_path + RAINfile)                   #name of mean pf precipitation rate at surface
rainf=np.squeeze(fninr['mtpr'][::])
rainf=np.transpose(rainf,(2,1,0))
rainf[np.where(np.abs(rainf)<1.e-7)] = 0



fnindl=dat(input_path + DLWRADfile)                  #downward long radiation flux
dlwrf=np.squeeze(fnindl["msdwlwrf"][::])
dlwrf=np.transpose(dlwrf,(2,1,0))


fninst=dat(input_path + TAIRfile)                  #surface air temperature
air=np.squeeze(fninst["t2m"][::])
air=np.transpose(air,(2,1,0))
air = air -273.15; #K to  C


fninsp=dat(input_path + PAIRfile)                  #surface air pressure
pres=np.squeeze(fninsp["sp"][::])
pres=np.transpose(pres,(2,1,0))
pres = pres*.01; #Pa to  mb



fninds=dat(input_path + DSWRADfile)                  #net shortwave radiation
swrad=np.squeeze(fninds["msnswrf"][::])
swrad=np.transpose(swrad,(2,1,0))


fninsh=dat(input_path + QAIRfile)                  #relative humidity (gfs)
shum=np.squeeze(fninsh["r"][::])
shum=np.transpose(shum,(2,1,0))


#####CONVERT SPECIFIC HUMIDITY TO RELATIVE HUMIDITY

# specific humidity at saturation (kg/kg).
# air temperature Ta (deg C). Dependence on air pressure, Pa, is small,

#    INPUT:   Ta - air temperature  [C]
#             Pa - (optional) pressure [mb]

#    OUTPUT:  q  - saturation specific humidity  [kg/kg]

#ew = 6.1121*(1.0007+3.46e-6*pres)*np.exp((17.502*air)/(240.97+air)); # in mb

#srs  = 0.62197*(ew/(pres-0.378*ew));                         # mb -> kg/kg

#shum=shum/srs
#shum=shum*100

###############################




tmf=fnin['time'][::]


# data de origem da simulacao

tzero=dates.datestr2num(str(Yorig)+'-'+str(1).zfill(2)+'-'+str(1).zfill(2) + ' ' +str(0).zfill(2)+':'+str(0).zfill(2)+':'+str(0).zfill(2))*24; #em horas                    

# passando para o msm referencial do matlab (0000,01,01)

tempo=tmf+ dates.datestr2num(str(1900)+'-'+str(1).zfill(2)+'-'+str(1).zfill(2))*24; #em horas | Reanalysis 2 #tempo do ncep no referencial do datenum

#data inicio da simulacao: (informacao que vem no romstools_param)

#time=(tempo-tzero)/24; #days since

time=(tempo-tzero)*(60*60); #seconds since ....


#acha o indice no qual a simulacao comeca;
first_time = time[0].copy()

i_dateref=np.where(time==first_time);


lnf=fnin['longitude'][::]
xf=np.where(lnf>180)
lnf[xf]=lnf[xf]-360; 
laf=fnin['latitude'][::]


[LNf,LAf]=np.meshgrid(lnf,laf);

LAf=np.flipud(LAf);

nt=(len(time)-(i_dateref[0][0]+ 1)); 

fa=4./24.; #freq de amostragem em dias (no caso da reanalise 6 por dia = 6/24)
na=nt*fa+1; ###num de dias para cycle_length !!!! Atencao !!!! 

L=len(lnf);
M=len(laf);

fname='u_era'+vari+'.nc';
print('    Cria ',fname)


type = 'BULK file from NCEP' ; 
history = 'ROMS' ;


g=ncgene(fname,'Uwind',L,M,nt,LAf,LNf,time, i_dateref, uwindf,'wind_time')    #u wind
g.start()


fname='v_era'+vari+'.nc';        #vwind
g.change_var(fname,'Vwind', vwindf, 'wind_time')
g.start()


fname='rain_era'+vari+'.nc';        #rain
g.change_var(fname,'rain', rainf, 'rain_time')
g.start()



fname='lwrad_era'+vari+'.nc';        #long wave radiation (using downward longwave radiation flux, not net)
g.change_var(fname,'lwrad_down', dlwrf, 'lrf_time')
g.start()


fname='temp_era'+vari+'.nc';        # surface air temperature
g.change_var(fname,'Tair', air, 'tair_time')
g.start()


fname='pressure_era'+vari+'.nc';        # surface air pressure
g.change_var(fname,'Pair', pres, 'pair_time')
g.start()



fname='srad_era'+vari+'.nc';        # net shortwave radiaton
g.change_var(fname,'swrad', swrad, 'srf_time')
g.start()


fname='humi_era'+vari+'.nc';        # specific humidity
g.change_var(fname,'Qair', shum, 'qair_time')
g.start()




