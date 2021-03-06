import pandas as pd
from matplotlib import dates
import numpy as np

def vel_conv(vel,dir):
  if dir <= 90:
    u = vel*np.sin(np.radians(dir))
    v = vel*np.cos(np.radians(dir))
  if dir > 90 and dir <=180:
    dir=dir-90
    u = vel*np.cos(np.radians(dir))
    v = -vel*np.sin(np.radians(dir))
  if dir > 180 and dir <=270:
    dir=dir-180
    u = -vel*np.sin(np.radians(dir))
    v = -vel*np.cos(np.radians(dir))
  if dir > 270 and dir <=360:
    dir=dir-270
    u = -vel*np.cos(np.radians(dir))
    v = vel*np.sin(np.radians(dir))
  return(u,v)  
  
###############adcp

adcphis=pd.read_csv('historico_itajai.txt')
#adcphis=pd.read_csv('historico_vitoria.txt')
#adcphis=pd.read_csv('historico_cabofrio.txt')
adcphis[adcphis.Lat<-100]=np.nan
adcphis = adcphis.dropna()
lat=adcphis['Lat'].min()
lon=adcphis['Lon'].min()

adcp=pd.read_csv('adcptratados_itajai.csv')
#adcp=pd.read_csv('adcptratados_vitoria_2.csv')
#adcp=pd.read_csv('adcptratados_cabofrio2.csv')
adcp['Lat']=lat
adcp['Lon']=lon
adcp['datas']=dates.datestr2num(adcp.data)

vel, dir='Cvel3','Cdir3'
#vel, dir='Cvel12','Cdir12'
#vel, dir='Cvel14','Cdir14'
#vel, dir='Cvel9','Cdir9'
#vel, dir='Cvel20','Cdir20'


vel=adcp[vel]
dir=adcp[dir]

u=np.zeros(adcp.shape[0])
v=np.zeros(adcp.shape[0])

for i in range(adcp.shape[0]):
  u[i],v[i] = vel_conv(vel[i], dir[i])
  u[i]=u[i]/1000
  v[i]=v[i]/1000

adcp['umed']=u
adcp['vmed']=v

#########################


from vgrid import s_coordinate, s_coordinate_2, s_coordinate_4
from netCDF4 import Dataset 
from scipy.interpolate import griddata, interp1d
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.pylab import *
from matplotlib import dates
from mpl_toolkits.axes_grid1 import make_axes_locatable
import glob
import xarray as xr

#avgfile=Dataset('HIS_FILE_20200421_5D0-20200428_5D0_hind_correct_year_WEAK_menor_azul_nopline_0005.nc') ##Vitoria 4,5
#avgfile=Dataset('HIS_FILE_rotate_cabofrio_3_0004.nc')        #CF 1/12
#avgfile=Dataset('HIS_FILE_rotate_cabofrio_3_HIGH_tmz_2_mm_22_0004.nc')               #CF1/36

lista = sorted(glob.glob('HIS_FILE_rotate_4_SUL_2_3_2_000[3-4]*')) #SUL 1/12

#lista = sorted(glob.glob('HIS_FILE_rotate_4_SUL_2_3_2_NEST_grid2_lp_all_small_ultra_smaler_000[3-4]*')) #SUL 1/36

avgfile = xr.open_mfdataset(lista, concat_dim='ocean_time')

hh=np.array(avgfile['ocean_time'][:].dt.strftime('%Y-%m-%d  %H:00:00'))

hhdates=dates.datestr2num(hh)
#fname_grd = 'azul_grd_era_NEW_menor_azul.nc'  #Vitoria

#fname_grd = 'rotate_cf_16_09.nc'    #CF 1/12
#fname_grd = 'CF_tmz_small026_2.nc'     #CF 1/36

fname_grd = 'grid_rotated_SUL_2.nc' #SUL 1/12
#fname_grd = 'grid_rotated_SUL_2_NEST_smaler.nc' #SUL 1/12


## Load ROMS grid.
grd = Dataset(fname_grd)
x_roms = grd.variables['lon_rho'][:]
y_roms = grd.variables['lat_rho'][:]
msk_roms = grd.variables['mask_rho'][:]
msk_romsv = grd.variables['mask_v'][:]
h_roms = grd.variables['h'][:]
sh2 = msk_roms.shape
etamax, ximax = sh2

ltlatlon=abs(y_roms - lat)
lglatlon=abs(x_roms - lon)
latlon=ltlatlon + lglatlon
ltmin=int(np.where(latlon == latlon.min())[0])
lgmin=int(np.where(latlon == latlon.min())[1])

theta_b = 0.4
theta_s = 5.0
tcline = 100.
klevels = 40
Vtransform = 2
Vstretching = 4
Spherical = True

if Vstretching==4:
  scoord = s_coordinate_4(h_roms, theta_b, theta_s, tcline, klevels)   #zeta is not used in the computation 
elif Vstretching==2:
  scoord = s_coordinate_2(h_roms, theta_b, theta_s, tcline, klevels)
elif Vstretching==1:
  scoord = s_coordinate(h_roms, theta_b, theta_s, tcline, klevels)

zr = -scoord.z_r[:]
 
zc=np.array([12.5])
#zc=np.array([47.5])
#zc=np.array([44])
#zc=np.array([33.5])
#zc=np.array([72.])


zc=zc[::-1]

uavg=np.array(avgfile['u_eastward'][:])

vavg=np.array(avgfile['v_northward'][:])

#time=avgfile['ocean_time'][:]/(24*60*60)

tempavg=np.array(avgfile['temp'][:])


ktl=3

uavg=uavg[:,:,ltmin-ktl:ltmin+ktl, lgmin-ktl:lgmin+ktl]
vavg=vavg[:,:,ltmin-ktl:ltmin+ktl, lgmin-ktl:lgmin+ktl]
tempavg=tempavg[:,:,ltmin-ktl:ltmin+ktl, lgmin-ktl:lgmin+ktl]
x_roms=x_roms[ltmin-ktl:ltmin+ktl, lgmin-ktl:lgmin+ktl]
y_roms=y_roms[ltmin-ktl:ltmin+ktl, lgmin-ktl:lgmin+ktl]

time=np.concatenate([time[:]/(24*60*60),time2[:]/(24*60*60)])

intu=np.zeros([uavg.shape[0],len(zc),x_roms.shape[0], x_roms.shape[1]])

intv=np.zeros([uavg.shape[0],len(zc),x_roms.shape[0], x_roms.shape[1]])

itemp=np.zeros([uavg.shape[0],len(zc),x_roms.shape[0], x_roms.shape[1]])

#zeta=avgfile['zeta'][:]

for j in range(intu.shape[2]):
  for k in range(intu.shape[3]):
    if (zr[-1,j,k] > zc.min()):
      zr[-1,j,k] = zc.min()

zr=zr[:,ltmin-ktl:ltmin+ktl, lgmin-ktl:lgmin+ktl]

UNDEF=np.nan

for i in range(intu.shape[0]):
  for j in range(intu.shape[2]):
    for k in range(intu.shape[3]):
      intu[i,:,j,k] = np.interp(-zc, -zr[:,j,k], uavg[i,:,j,k], right=UNDEF, left=UNDEF)
      intv[i,:,j,k] = np.interp(-zc, -zr[:,j,k], vavg[i,:,j,k], right=UNDEF, left=UNDEF)
#      itemp[i,:,j,k] = np.interp(-zc, -zr[:,j,k], tempavg[i,:,j,k], right=UNDEF, left=UNDEF)

intu=np.squeeze(intu)
intv=np.squeeze(intv)
itemp=np.squeeze(itemp)


#begindate=avgfile['ocean_time'].units[14:]
#begindate=dates.datestr2num(begindate)
romstime=hhdates

########adcp

ff=np.where(  (adcp['datas'] > romstime.min()) & (adcp['datas'] < romstime.max()) )

newadcp=adcp.loc[ff]

timeostia=newadcp['datas']

############################################################

intunewtime=np.zeros([len(timeostia),x_roms.shape[0], x_roms.shape[1]])

intvnewtime=np.zeros([len(timeostia),x_roms.shape[0], x_roms.shape[1]])

tempnewtime=np.zeros([len(timeostia),x_roms.shape[0], x_roms.shape[1]])


for j in range(intu.shape[1]):
  for k in range(intu.shape[2]):
    intunewtime[:,j,k] = np.interp(timeostia, romstime, intu[:,j,k], right=UNDEF, left=UNDEF)
    intvnewtime[:,j,k] = np.interp(timeostia, romstime, intv[:,j,k], right=UNDEF, left=UNDEF)
    tempnewtime[:,j,k] = np.interp(timeostia, romstime, itemp[:,j,k], right=UNDEF, left=UNDEF)

intunewtime[intunewtime>100]=np.nan
#plt.pcolor(x_roms, y_roms,np.ma.masked_invalid(intunewtime[-1,:]));plt.show()

vsitu=np.zeros([len(timeostia)])
usitu=np.zeros([len(timeostia)])

for i in range(len(vsitu)):
  vsitu[i]=griddata((x_roms.ravel(),y_roms.ravel()), intvnewtime[i,:].ravel(), (lon,lat))
  usitu[i]=griddata((x_roms.ravel(),y_roms.ravel()), intunewtime[i,:].ravel(), (lon,lat))

vmed=newadcp['vmed']

umed=newadcp['umed']

timevec=newadcp['datas'].values

vbeca=vsitu.copy()
ubeca=usitu.copy()

vson=vsitu.copy()
uson=usitu.copy()



vbeca=vbeca[0:1415]
ubeca=ubeca[0:1415]

#######################################################################  MERCATOR

from vgrid import s_coordinate, s_coordinate_2, s_coordinate_4
from netCDF4 import Dataset 
from scipy.interpolate import griddata, interp1d
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.pylab import *
from matplotlib import dates
from mpl_toolkits.axes_grid1 import make_axes_locatable
       

#avgfile=Dataset('HIS_FILE_20200421_5D0-20200428_5D0_hind_correct_year_WEAK_menor_azul_nopline_0005.nc') ##Vitoria 4,5
#avgfile=Dataset('HIS_FILE_rotate_cabofrio_3_2_0004.nc')   #cabo frio 3,4


lista = sorted(glob.glob('HIS_FILE_rotate_4_SUL_2_3_2_NEST_grid2_lp_all_small_ultra_smaler_000[3-4]*'))

avgfile = xr.open_mfdataset(lista, concat_dim='ocean_time')

hh=np.array(avgfile['ocean_time'][:].dt.strftime('%Y-%m-%d  %H:00:00'))

hhdates=dates.datestr2num(hh)

#fname_grd = 'azul_grd_era_NEW_menor_azul.nc'  #Vitoria
#fname_grd = 'rotate_cf_16_09_2.nc'    #CF
fname_grd = 'grid_rotated_SUL_2_NEST_smaler.nc' #SUL 1/36


#time=avgfile['ocean_time'][:]/(24*60*60)
#begindate=avgfile['ocean_time'].units[14:]
#begindate=dates.datestr2num(begindate)
#romstime=begindate+time
romstime=hhdates

timeini=dates.num2date(romstime[0]  - 1).strftime("%Y%m%d")
timeend=dates.num2date(romstime[-1]).strftime("%Y%m%d")

input_path='/home/fernando/roms/src/Projects/hindcast_2/mercator/MYOCEAN_AZUL_FORECAST_'    #
#input_path='/home/fernando/roms/src/Projects/hindcast/mercator/MYOCEAN_AZUL_FORECAST_'    #


a=[timeini, timeend]

numdays= np.diff(dates.datestr2num(a)) + 2

n_time=int(numdays[0])

a = dates.datestr2num(a)

time_vec = np.arange(numdays) + a[0]

romstimer=romstime.copy()

romstime=time_vec.copy()


rt=0
file=Dataset(input_path + dates.num2date( time_vec[rt]).strftime("%Y%m%d")+'.nc')
ut=file['uo'][:]

uavg=np.zeros([int(numdays)] + list(np.squeeze(ut).shape))
vavg=np.zeros([int(numdays)] + list(np.squeeze(ut).shape))


for rt in range(int(numdays)):
  #load HYCOM DATA
  nc_hycom=input_path + dates.num2date( time_vec[rt]).strftime("%Y%m%d")+'.nc'
  file=Dataset(nc_hycom)
  uavg[rt,:]=np.squeeze(file['uo'][:])
  vavg[rt,:]=np.squeeze(file['vo'][:])

uavg[uavg<-100]=np.nan
vavg[vavg<-100]=np.nan
zr = file['depth'][:]

zc=np.array([12.5])

#zc=np.array([47.5])

#zc=np.array([51.0])

#zc=np.array([44])

#zc=np.array([72.])



x_fm=file['longitude'][:]
y_fm=file['latitude'][:]
x_fm,y_fm=np.meshgrid(x_fm,y_fm)

minlon = x_fm[0,:] - lon
iml = np.where(np.absolute(minlon)==np.absolute(minlon).min())[0][0]
maxlon = x_fm[0,:] - lon
imxl = np.where(np.absolute(maxlon)==np.absolute(maxlon).min())[0][0]

minlat = y_fm[:,0] - lat
imla = np.where(np.absolute(minlat)==np.absolute(minlat).min())[0][0]
maxlat = y_fm[:,0] - lat
imxla = np.where(np.absolute(maxlat)==np.absolute(maxlat).min())[0][0]

lim=3

uavg=uavg[:,:,imxla-lim:imla+lim,iml-lim:imxl+lim]
vavg=vavg[:,:,imxla-lim:imla+lim,iml-lim:imxl+lim]
tempavg=tempavg[:,:imxla-lim:imla+lim,iml-lim:imxl+lim]


x_fm = x_fm[imxla-lim:imla+lim,iml-lim:imxl+lim]

y_fm = y_fm[imxla-lim:imla+lim,iml-lim:imxl+lim]


intu=np.zeros([uavg.shape[0],len(zc),uavg.shape[2], uavg.shape[3]])

intv=np.zeros([uavg.shape[0],len(zc),uavg.shape[2], uavg.shape[3]])

itemp=np.zeros([uavg.shape[0],len(zc),uavg.shape[2], uavg.shape[3]])

UNDEF=np.nan

for i in range(intu.shape[0]):
  for j in range(intu.shape[2]):
    for k in range(intu.shape[3]):
      intu[i,:,j,k] = np.interp(zc, zr, uavg[i,:,j,k], right=UNDEF, left=UNDEF)
      intv[i,:,j,k] = np.interp(zc, zr, vavg[i,:,j,k], right=UNDEF, left=UNDEF)
#      itemp[i,:,j,k] = np.interp(zc, zr, tempavg[i,:,j,k], right=UNDEF, left=UNDEF)

intu=np.squeeze(intu)
intv=np.squeeze(intv)
itemp=np.squeeze(itemp)


import pandas as pd
from matplotlib import dates
import numpy as np

def vel_conv(vel,dir):
  if dir <= 90:
    u = vel*np.sin(np.radians(dir))
    v = vel*np.cos(np.radians(dir))
  if dir > 90 and dir <=180:
    dir=dir-90
    u = vel*np.cos(np.radians(dir))
    v = -vel*np.sin(np.radians(dir))
  if dir > 180 and dir <=270:
    dir=dir-180
    u = -vel*np.sin(np.radians(dir))
    v = -vel*np.cos(np.radians(dir))
  if dir > 270 and dir <=360:
    dir=dir-270
    u = -vel*np.cos(np.radians(dir))
    v = vel*np.sin(np.radians(dir))
  return(u,v)  
  
###############adcp

adcphis=pd.read_csv('historico_itajai.txt')
#adcphis=pd.read_csv('historico_vitoria.txt')
#adcphis=pd.read_csv('historico_cabofrio.txt')
adcphis[adcphis.Lat<-100]=np.nan
adcphis = adcphis.dropna()
lat=adcphis['Lat'].min()
lon=adcphis['Lon'].min()

adcp=pd.read_csv('adcptratados_itajai.csv')
#adcp=pd.read_csv('adcptratados_vitoria_2.csv')
#adcp=pd.read_csv('adcptratados_cabofrio2.csv')
adcp['Lat']=lat
adcp['Lon']=lon
adcp['datas']=dates.datestr2num(adcp.data)

vel, dir='Cvel3','Cdir3'
#vel, dir='Cvel12','Cdir12'
#vel, dir='Cvel14','Cdir14'
#vel, dir='Cvel9','Cdir9'
#vel, dir='Cvel20','Cdir20'

vel=adcp[vel]
dir=adcp[dir]

u=np.zeros(adcp.shape[0])
v=np.zeros(adcp.shape[0])

for i in range(adcp.shape[0]):
  u[i],v[i] = vel_conv(vel[i], dir[i])
  u[i]=u[i]/1000
  v[i]=v[i]/1000

adcp['umed']=u
adcp['vmed']=v


ff=np.where(  (adcp['datas'] > romstimer.min()) & (adcp['datas'] < romstimer.max()) )

newadcp=adcp.loc[ff]

timeostia=newadcp['datas']

############################################################

intunewtimer=np.zeros([len(romstimer),intu.shape[1], intu.shape[2]])

intvnewtimer=np.zeros([len(romstimer),intu.shape[1], intu.shape[2]])

itempnewtimer=np.zeros([len(romstimer),intu.shape[1], intu.shape[2]])


for j in range(intu.shape[1]):
  for k in range(intu.shape[2]):
    intunewtimer[:,j,k] = np.interp(romstimer, romstime, intu[:,j,k], right=UNDEF, left=UNDEF)
    intvnewtimer[:,j,k] = np.interp(romstimer, romstime, intv[:,j,k], right=UNDEF, left=UNDEF)
    itempnewtimer[:,j,k] = np.interp(romstimer, romstime, itemp[:,j,k], right=UNDEF, left=UNDEF)


intunewtime=np.zeros([len(timeostia),intu.shape[1], intu.shape[2]])

intvnewtime=np.zeros([len(timeostia),intu.shape[1], intu.shape[2]])

tempnewtime=np.zeros([len(timeostia),intu.shape[1], intu.shape[2]])


for j in range(intu.shape[1]):
  for k in range(intu.shape[2]):
    intunewtime[:,j,k] = np.interp(timeostia, romstimer, intunewtimer[:,j,k], right=UNDEF, left=UNDEF)
    intvnewtime[:,j,k] = np.interp(timeostia, romstimer, intvnewtimer[:,j,k], right=UNDEF, left=UNDEF)
    tempnewtime[:,j,k] = np.interp(timeostia, romstimer, itempnewtimer[:,j,k], right=UNDEF, left=UNDEF)

intunewtime[intunewtime>100]=np.nan

#plt.pcolor(np.ma.masked_invalid(intunewtime[-1,:]));plt.show()

vsitu=np.zeros([len(timeostia)])
usitu=np.zeros([len(timeostia)])


for i in range(len(vsitu)):
  vsitu[i]=griddata((x_fm.ravel(),y_fm.ravel()), intvnewtime[i,:].ravel(), (lon,lat))
  usitu[i]=griddata((x_fm.ravel(),y_fm.ravel()), intunewtime[i,:].ravel(), (lon,lat))


vmed=newadcp['vmed']

umed=newadcp['umed']

timevec=newadcp['datas'].values

from utils import weim





vmedt=weim(vmed,81)
umedt=weim(umed,81)
vsitut=weim(vsitu,31)
usitut=weim(usitu,31)
vsitutdad=weim(vbeca,61)
usitutdad=weim(ubeca,61)
vsitutnest=weim(vson, 61)
usitutnest=weim(uson, 61)


umedt=weim(umed,81)
usitut=weim(usitu,81)



vmedt=weim(vmed1,91)
vsitut=weim(vsitu1,61)


umedt=weim(umed1,91)
usitut=weim(usitu1,61)



plt.plot(dates.num2date(timevec),vmedt, 'red', label='PNBOIA')
plt.plot(dates.num2date(timevec),vsitut, 'blue', label='ROMS')
ax = plt.gca()
ax.legend(loc='best')
plt.show()


plt.plot(dates.num2date(timevec),umedt, 'red',label='PNBOIA')
plt.plot(dates.num2date(timevec),usitut, 'blue',  label='ROMS')
ax = plt.gca()
ax.legend(loc='best')
plt.show()

from scipy import stats
from sklearn.metrics import mean_squared_error as mse

stats.pearsonr(vmedt,vsitut)

rmsev=np.sqrt(mse(vmedt,vsitut))
print(rmsev)


stats.pearsonr(umedt,usitut)

rmseu=np.sqrt(mse(umedt,usitut))
print(rmseu)



from scipy import stats
from sklearn.metrics import mean_squared_error as mse

stats.pearsonr(vmedt,vsitutnest)

rmsev=np.sqrt(mse(vmedt,vsitutnest))
print(rmsev)


stats.pearsonr(umedt,usitutnest)

rmseu=np.sqrt(mse(umedt,usitutnest))
print(rmseu)


from scipy import stats
from sklearn.metrics import mean_squared_error as mse

stats.pearsonr(vmedt,vsitutdad)

rmsev=np.sqrt(mse(vmedt,vsitutdad))
print(rmsev)


stats.pearsonr(umedt,usitutdad)

rmseu=np.sqrt(mse(umedt,usitutdad))
print(rmseu)


vmedt=weim(vmed,81)
umedt=weim(umed,81)
vsitut=weim(vsitu,31)
usitut=weim(usitu,31)
vsitutdad=weim(vbeca,61)
usitutdad=weim(ubeca,61)
vsitutnest=weim(vson, 61)
usitutnest=weim(uson, 61)

plt.style.use('ggplot')

fig, axs = plt.subplots(2, sharex=True, figsize=(9,6))
#fig.suptitle('Vertically stacked subplots')
axs[0].plot(dates.num2date(timevec),vsitut, 'r--', label='V-component MERCATOR')
axs[0].plot(dates.num2date(timevec),vsitutdad, 'green',linestyle='--', label='V-component ROMS PARENT')
#axs[0].plot(dates.num2date(timevec),vsitutdad, 'green',linestyle='--', label='V-component ROMS')
axs[0].plot(dates.num2date(timevec),vsitutnest, 'black',linestyle='--', label='V-component ROMS NEST')
axs[0].plot(dates.num2date(timevec),vmedt, 'blue', label='V-component PNBOIA')
legend=axs[0].legend(loc=2, fontsize='xx-small')
legend.get_frame().set_facecolor('grey')
axs[0].set_ylim([-0.3, 0.3])
axs[0].tick_params(labelsize=8)



axs[1].plot(dates.num2date(timevec),usitut, 'r--',  label='U-component MERCATOR')
axs[1].plot(dates.num2date(timevec),usitutdad, 'green',linestyle='--',  label='U-component ROMS PARENT')
#axs[1].plot(dates.num2date(timevec),usitutdad, 'green',linestyle='--',  label='U-component ROMS')
axs[1].plot(dates.num2date(timevec),usitutnest, 'black',linestyle='--', label='U-component ROMS NEST')
axs[1].plot(dates.num2date(timevec),umedt, 'blue',label='U-component PNBOIA')
axs[1].legend(loc=2)
legend=axs[1].legend(loc=2, fontsize='xx-small')
legend.get_frame().set_facecolor('grey')
axs[1].set_ylim([-0.5, 0.3])
axs[1].tick_params(labelsize=8)


fig.text(0.06, 0.5, 'Velocity (m/s)', ha='center', va='center', rotation='vertical', fontsize=11)
plt.xticks(fontsize=7,rotation=35)
plt.gcf().subplots_adjust(bottom=0.15)
#fig.text(0.5, 0.9,'Depth = ' + str(float(zc)) + ' m',  fontsize=12, ha='center', va='center')
fig.text(0.5, 0.9,'Depth = ' + '5.5' + ' m',  fontsize=12, ha='center', va='center')

plt.show()






fig, axs = plt.subplots(2, sharex=True, figsize=(8,5))
#fig.suptitle('Vertically stacked subplots')
axs[0].plot(dates.num2date(timevec),vsitut, 'red', label='V-component MERCATOR')
axs[0].plot(dates.num2date(timevec),vmedt, 'blue', label='V-component PNBOIA')
legend=axs[0].legend(loc=2, fontsize='x-small')
legend.get_frame().set_facecolor('grey')
axs[0].set_ylim([-0.3, 0.3])


axs[1].plot(dates.num2date(timevec),usitut, 'red',  label='U-componen MERCATOR')
axs[1].plot(dates.num2date(timevec),umedt, 'blue',label='U-component PNBOIA')
axs[1].legend(loc=2)
legend=axs[1].legend(loc=2, fontsize='x-small')
legend.get_frame().set_facecolor('grey')
axs[1].set_ylim([-0.2, 0.2])
fig.text(0.06, 0.5, 'Velocity (m/s)', ha='center', va='center', rotation='vertical', fontsize=12)


plt.xticks(fontsize=8,rotation=35)
plt.gcf().subplots_adjust(bottom=0.15)

plt.show()


