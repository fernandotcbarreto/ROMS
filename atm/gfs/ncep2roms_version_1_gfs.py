from parameters import *
from datetime import date
from netCDF4 import Dataset as dat
import numpy as np
from forclass import ncgene
import matplotlib.dates as dates
import matplotlib.pyplot as plt

UWINDfile='GFS_F_UWIND_20200415.nc'
VWINDfile='GFS_F_VWIND_20200415.nc'
RAINfile='GFS_F_RAIN_20200415.nc'
DLWRADfile='GFS_F_DLWRAD_20200415.nc'
TAIRfile='GFS_F_TAIR_20200415.nc'
PAIRfile='GFS_F_PAIR_20200415.nc'
DSWRADfile='GFS_F_DSWRAD_20200415.nc'
USWRADfile='GFS_F_USWRAD_20200415.nc'
QAIRfile='GFS_F_QAIR_20200415.nc'

fnin=dat(input_path + UWINDfile)                  #name of U wind file
uwindf=np.squeeze(fnin["ugrd10m"][::])
uwindf=np.transpose(uwindf,(2,1,0))


fninv=dat(input_path + VWINDfile)                  #name of V wind file
vwindf=np.squeeze(fninv["vgrd10m"][::])
vwindf=np.transpose(vwindf,(2,1,0))


fninr=dat(input_path + RAINfile)                   #name of mean pf precipitation rate at surface
rainf=np.squeeze(fninr['pratesfc'][::])
rainf=np.transpose(rainf,(2,1,0))
rainf[np.where(np.abs(rainf)<1.e-7)] = 0



fnindl=dat(input_path + DLWRADfile)                  #downward long radiation flux
dlwrf=np.squeeze(fnindl["dlwrfsfc"][::])
dlwrf=np.transpose(dlwrf,(2,1,0))


fninst=dat(input_path + TAIRfile)                  #surface air temperature
air=np.squeeze(fninst["tmp2m"][::])
air=np.transpose(air,(2,1,0))
air = air -273.15; #K to  C


fninsp=dat(input_path + PAIRfile)                  #surface air pressure
pres=np.squeeze(fninsp["pressfc"][::])
pres=np.transpose(pres,(2,1,0))
pres = pres*.01; #Pa to  mb



fninds=dat(input_path + DSWRADfile)                  #downward shortwave radiation
dshor=np.squeeze(fninds["dswrfsfc"][::])
dshor=np.transpose(dshor,(2,1,0))

fninus=dat(input_path + USWRADfile)                  #upward shortwave radiation
ushor=np.squeeze(fninus["uswrfsfc"][::])
ushor=np.transpose(ushor,(2,1,0))

swrad=dshor-ushor;                          #net shortwave radiation



fninsh=dat(input_path + QAIRfile)                  #relative humidity (gfs)
shum=np.squeeze(fninsh["rh2m"][::])
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


timeref = '2013-01-01 00:00:00'                   # Reference time at bry file. Azul project is seconds since 2013-01-01 00:00:00]


tmf=fnin['time'][::]

yearf = int(fnin['time'].units[13:18])
monf  = int(fnin['time'].units[19:21])
dayf  = int(fnin['time'].units[22:24])

Ymin = yearf
Mmin = monf
Dmin = dayf
# data de origem da simulacao

tzero=dates.datestr2num(str(Yorig)+'-'+str(1).zfill(2)+'-'+str(1).zfill(2) + ' ' +str(0).zfill(2)+':'+str(0).zfill(2)+':'+str(0).zfill(2))*24; #em horas                    

# passando para o msm referencial do matlab (0000,01,01)

first_time = dates.datestr2num(str(Ymin)+'-'+str(Mmin).zfill(2)+'-'+str(Dmin).zfill(2) + ' ' +str(Hmin).zfill(2)+':'+str(0).zfill(2)+':'+str(0).zfill(2)) -tzero/24;    #convert to julian day

tempo=tmf/60.+ dates.datestr2num(str(yearf)+'-'+str(monf).zfill(2)+'-'+str(dayf).zfill(2))*24; #em horas | Reanalysis 2 #tempo do ncep no referencial do datenum

#data inicio da simulacao: (informacao que vem no romstools_param)

#time=(tempo-tzero)/24; days since
time=(tempo-tzero)*(60*60); #seconds since ....


#vetor para o roms, deve ser 0 no instante de inicio da simulacao. (definido no
#Dateref = 0)

#acha o indice no qual a simulacao comeca;
first_time = first_time*(24*60*60)


i_dateref=np.where(time==first_time);


lnf=fnin['lon'][::]
xf=np.where(lnf>180)
lnf[xf]=lnf[xf]-360; 
laf=fnin['lat'][::]


[LNf,LAf]=np.meshgrid(lnf,laf);

#LAf=np.flipud(LAf);   #GFS already comes in right lat orientation (south as 0 index)

nt=(len(time)-(i_dateref[0][0]+ 1)); 

fa=4./24.; #freq de amostragem em dias (no caso da reanalise 6 por dia = 6/24)
na=nt*fa+1; ###num de dias para cycle_length !!!! Atencao !!!! 

L=len(lnf);
M=len(laf);

fname='20131_u_gfs.nc';
print('    Cria ',fname)


type = 'BULK file from NCEP' ; 
history = 'ROMS' ;


g=ncgene(fname,'Uwind',L,M,nt,LAf,LNf,time, i_dateref, uwindf,'wind_time')    #u wind
g.start()


fname='20131_v_gfs.nc';        #vwind
g.change_var(fname,'Vwind', vwindf, 'wind_time')
g.start()


fname='20131_rain_gfs.nc';        #rain
g.change_var(fname,'rain', rainf, 'rain_time')
g.start()



fname='20131_lwrad_gfs.nc';        #long wave radiation (using downward longwave radiation flux, not net)
g.change_var(fname,'lwrad_down', dlwrf, 'lrf_time')
g.start()


fname='20131_temp_gfs.nc';        # surface air temperature
g.change_var(fname,'Tair', air, 'tair_time')
g.start()


fname='20131_pressure_gfs.nc';        # surface air pressure
g.change_var(fname,'Pair', pres, 'pair_time')
g.start()



fname='20131_srad_gfs.nc';        # net shortwave radiaton
g.change_var(fname,'swrad', swrad, 'srf_time')
g.start()


fname='20131_humi_corrigido_gfs.nc';        # specific humidity
g.change_var(fname,'Qair', shum, 'qair_time')
g.start()