from parameters import *
from datetime import date
from netCDF4 import Dataset as dat
import numpy as np
from forclass import ncgene
import matplotlib.dates as dates
import matplotlib.pyplot as plt

fnin=dat(input_path + 'u_2013.nc')                  #name of U wind file
uwindf=np.squeeze(fnin["uwnd"][::])
uwindf=np.transpose(uwindf,(2,1,0))


fninv=dat(input_path + 'v_2013.nc')                  #name of V wind file
vwindf=np.squeeze(fninv["vwnd"][::])
vwindf=np.transpose(vwindf,(2,1,0))


fninr=dat(input_path+'rain_2013.nc')                   #name of mean pf precipitation rate at surface
rainf=np.squeeze(fninr['prate'][::])
rainf=np.transpose(rainf,(2,1,0))
rainf[np.where(np.abs(rainf)<1.e-7)] = 0



fnindl=dat(input_path+'dlrf_2013.nc')                  #downward long radiation flux
dlwrf=np.squeeze(fnindl["dlwrf"][::])
dlwrf=np.transpose(dlwrf,(2,1,0))


fninst=dat(input_path+'temp_2013.nc')                  #surface air temperature
air=np.squeeze(fninst["air"][::])
air=np.transpose(air,(2,1,0))
air = air -273.15; #K to  C


fninsp=dat(input_path+'pres_2013.nc')                  #surface air pressure
pres=np.squeeze(fninsp["pres"][::])
pres=np.transpose(pres,(2,1,0))
pres = pres*.01; #Pa to  mb



fninds=dat(input_path + 'downshort_2013.nc')                  #downward shortwave radiation
dshor=np.squeeze(fninds["dswrf"][::])
dshor=np.transpose(dshor,(2,1,0))

fninus=dat(input_path + 'upshor_2013.nc')                  #upward shortwave radiation
ushor=np.squeeze(fninus["uswrf"][::])
ushor=np.transpose(ushor,(2,1,0))

swrad=dshor-ushor;                          #net shortwave radiation



fninsh=dat(input_path+'humi_2013.nc')                  #specific humidity
shum=np.squeeze(fninsh["shum"][::])
shum=np.transpose(shum,(2,1,0))


#####CONVERT SPECIFIC HUMIDITY TO RELATIVE HUMIDITY

# specific humidity at saturation (kg/kg).
# air temperature Ta (deg C). Dependence on air pressure, Pa, is small,

#    INPUT:   Ta - air temperature  [C]
#             Pa - (optional) pressure [mb]

#    OUTPUT:  q  - saturation specific humidity  [kg/kg]

ew = 6.1121*(1.0007+3.46e-6*pres)*np.exp((17.502*air)/(240.97+air)); # in mb

srs  = 0.62197*(ew/(pres-0.378*ew));                         # mb -> kg/kg

shum=shum/srs
shum=shum*100

###############################




tmf=fnin['time'][::]


# data de origem da simulacao

tzero=dates.datestr2num(str(Yorig)+'-'+str(1).zfill(2)+'-'+str(1).zfill(2) + ' ' +str(0).zfill(2)+':'+str(0).zfill(2)+':'+str(0).zfill(2))*24; #em horas                    

# passando para o msm referencial do matlab (0000,01,01)

tempo=tmf+ dates.datestr2num(str(1800)+'-'+str(1).zfill(2)+'-'+str(1).zfill(2))*24; #em horas | Reanalysis 2 #tempo do ncep no referencial do datenum

#data inicio da simulacao: (informacao que vem no romstools_param)

#time=(tempo-tzero)/24; #days since

time=(tempo-tzero)*(60*60); #seconds since ....


#acha o indice no qual a simulacao comeca;
first_time = time[0].copy()

i_dateref=np.where(time==first_time);


lnf=fnin['lon'][::]
xf=np.where(lnf>180)
lnf[xf]=lnf[xf]-360; 
laf=fnin['lat'][::]


[LNf,LAf]=np.meshgrid(lnf,laf);

LAf=np.flipud(LAf);

nt=(len(time)-(i_dateref[0][0]+ 1)); 

fa=4./24.; #freq de amostragem em dias (no caso da reanalise 6 por dia = 6/24)
na=nt*fa+1; ###num de dias para cycle_length !!!! Atencao !!!! 

L=len(lnf);
M=len(laf);

fname='20131_u_11m.nc';
print('    Cria ',fname)


type = 'BULK file from NCEP' ; 
history = 'ROMS' ;


g=ncgene(fname,'Uwind',L,M,nt,LAf,LNf,time, i_dateref, uwindf,'wind_time')    #u wind
g.start()


fname='20131_v_11m.nc';        #vwind
g.change_var(fname,'Vwind', vwindf, 'wind_time')
g.start()


fname='20131_rain_11m.nc';        #rain
g.change_var(fname,'rain', rainf, 'rain_time')
g.start()



fname='20131_lwrad_11m.nc';        #long wave radiation (using downward longwave radiation flux, not net)
g.change_var(fname,'lwrad_down', dlwrf, 'lrf_time')
g.start()


fname='20131_temp_11m.nc';        # surface air temperature
g.change_var(fname,'Tair', air, 'tair_time')
g.start()


fname='20131_pressure_11m.nc';        # surface air pressure
g.change_var(fname,'Pair', pres, 'pair_time')
g.start()



fname='20131_srad_11m.nc';        # net shortwave radiaton
g.change_var(fname,'swrad', swrad, 'srf_time')
g.start()


fname='20131_humi_11m.nc';        # specific humidity
g.change_var(fname,'Qair', shum, 'qair_time')
g.start()




