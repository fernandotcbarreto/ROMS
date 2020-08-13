import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as dat
from mpl_toolkits.basemap import Basemap


file=dat('ocean_BRSE_avg_clim_nud_kate_gls_sem_spline_0001.nc')

time=file['ocean_time'][::]

temp=file['temp'][::]

lat=file['lat_rho'][::]

lon=file['lon_rho'][::]

pathfig='/home/fernando/roms/src/Projects/Azul/figs/'


for i in range(len(time)):
  fig, ax1 = plt.subplots()
  a=ax1.pcolor(temp[i,0,:,:], vmin=temp.min(), vmax=temp.max())
  cbar=fig.colorbar(a, shrink=0.9, extend='both')
  plt.savefig(pathfig+str(i)+'fundo.png')
  
  
 ############################################   2d maps
v=file['v_northward'][:]
u=file['u_eastward'][:]

tim=5

lay=0

sp=3

val=np.sqrt((u[tim,lay,:,:]**2)+(v[tim,lay,:,:]**2))

#a=plt.pcolor(lon,lat, np.squeeze(val),vmin=val.min(),vmax=val.max()) 
a=plt.pcolor(lon,lat, np.squeeze(val),vmin=val.min(),vmax=val.max())
plt.quiver(lon[0:-1:sp, 0:-1:sp],lat[0:-1:sp, 0:-1:sp],u[tim,lay,0:-1:sp, 0:-1:sp], v[tim,lay,0:-1:sp, 0:-1:sp])
plt.colorbar(a)
plt.show()



a=temp[-1,-1,:]
plt.pcolor(lon,lat, np.squeeze(a),vmin=a.min(),vmax=a.max())
plt.colorbar()
plt.show()

sp=2
tim=-1
lay=-1

rj=(-42.7,-22.5)

nlim=-18
slim=-26
wlim=-46
elim=-37

val=np.sqrt((u[tim,lay,:,:]**2)+(v[tim,lay,:,:]**2))


f, axs = plt.subplots(figsize=(10,10))

# set up orthographic map projection with
# perspective of satellite looking down at 50N, 100W.
# use low resolution coastlines.
map = Basemap(projection='merc', llcrnrlat=slim, urcrnrlat=nlim,llcrnrlon=wlim, urcrnrlon=elim,resolution='l')
# draw coastlines, country boundaries, fill continents.
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.fillcontinents(color='#ddaa66',lake_color='aqua')
parallels = np.arange(-30,-10,2)
map.drawparallels(parallels,labels=[1,0,0,1], linewidth=0.0)
meridians = np.arange(-55,-30,2)
map.drawmeridians(meridians,labels=[1,0,0,1], linewidth=0.0)
map.readshapefile('/mnt/c/Users/Fernando/Desktop/shape_unid_fed/lim_unidade_federacao_a', 'lim_unidade_federacao_a')
x, y = map(lon,lat)
a=axs.pcolor(x,y, np.squeeze(val),vmin=val.min(),vmax=val.max())
axs.quiver(x[0:-1:sp, 0:-1:sp],y[0:-1:sp, 0:-1:sp],u[tim,lay,0:-1:sp, 0:-1:sp], v[tim,lay,0:-1:sp, 0:-1:sp])
axs.text(map(rj[0], rj[1])[0],map(rj[0], rj[1])[1],'RJ', fontsize=14, fontweight='bold',
        bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 5})
plt.colorbar(a)
plt.show()



  
file=dat('ocean_BRSE_avg_0001.nc')

temp=file['zeta'][::]

time=file['ocean_time'][::]

lat=file['lat_rho'][::]

lon=file['lon_rho'][::]

pathfig='/home/fernando/roms/src/Projects/Azul/figs/'


for i in range(len(time)):
  fig, ax1 = plt.subplots()
  a=ax1.pcolor(temp[i,:,:], vmin=temp.min(), vmax=temp.max())
  cbar=fig.colorbar(a, shrink=0.9, extend='both')
  plt.savefig(pathfig+str(i)+'fundo.png')
  print(i)
  

  
  

  
  

azul_grd2.nc 

prooceano_myocean_mais_13_03_ini.nc

prooceano_myocean_cortado_1_12_bry.nc


azul_tides_otps_13_03.nc



file=dat('ocean_BRSE_avg_0001.nc')

temp=file['temp'][:]

lon=file['lon_rho'][:]

lat=file['lat_rho'][:]

temp[temp>50]=np.nan
temp[temp<-50]=np.nan

plt.pcolor(lon,lat,np.ma.masked_invalid(temp[1,0,::]));plt.colorbar();plt.show()

temp=file['u'][:]

lon=file['lon_rho'][:]

lat=file['lat_rho'][:]

temp[temp>50]=np.nan
temp[temp<-50]=np.nan

plt.pcolor(lon,lat,np.ma.masked_invalid(temp[1,0,::]));plt.colorbar();plt.show()


temp=file['salt'][:]

lon=file['lon_rho'][:]

lat=file['lat_rho'][:]

temp[temp>50]=np.nan
temp[temp<-50]=np.nan

plt.pcolor(lon,lat,np.ma.masked_invalid(temp[1,10,::]));plt.colorbar();plt.show()

temp=file['v'][:]

lon=file['lon_rho'][:]

lat=file['lat_rho'][:]

temp[temp>50]=np.nan
temp[temp<-50]=np.nan

plt.pcolor(lon,lat,np.ma.masked_invalid(temp[1,-1,::]));plt.colorbar();plt.show()


temp=file['zeta'][:]

lon=file['lon_rho'][:]

lat=file['lat_rho'][:]

temp[temp>50]=np.nan
temp[temp<-50]=np.nan

plt.pcolor(lon,lat,np.ma.masked_invalid(temp[1,::]));plt.colorbar();plt.show()



temp=file['ubar'][:]

lon=file['lon_rho'][:]

lat=file['lat_rho'][:]

temp[temp>50]=np.nan
temp[temp<-50]=np.nan

plt.pcolor(lon,lat,np.ma.masked_invalid(temp[1,::]));plt.colorbar();plt.show()

temp=file['vbar'][:]

lon=file['lon_rho'][:]

lat=file['lat_rho'][:]

temp[temp>50]=np.nan
temp[temp<-50]=np.nan

plt.pcolor(lon,lat,np.ma.masked_invalid(temp[1,::]));plt.colorbar();plt.show()



