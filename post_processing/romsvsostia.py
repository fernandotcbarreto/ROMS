#roms vs ostia

from vgrid import s_coordinate, s_coordinate_2, s_coordinate_4
from netCDF4 import Dataset 
from scipy.interpolate import griddata, interp1d
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.pylab import *
from matplotlib import dates
from mpl_toolkits.axes_grid1 import make_axes_locatable

#avgfile=Dataset('HIS_FILE_20201109_5D0-20201118_5D0_hind.nc')
avgfile=Dataset('his_no_ass.nc')
#avgfile=Dataset('BRSE_fwd_outer1.nc')
#avgfile2=Dataset('HIS_FILE_20200421_5D0-20200428_5D0_hind_correct_year_WEAK_menor_azul_nopline_0010.nc')

fname_grd = 'azul_grd2.nc'  ## Smooth topography grid.   


## Load ROMS grid.
grd = Dataset(fname_grd)
x_roms = grd.variables['lon_rho'][:]
y_roms = grd.variables['lat_rho'][:]
msk_roms = grd.variables['mask_rho'][:]
msk_romsv = grd.variables['mask_v'][:]
h_roms = grd.variables['h'][:]
sh2 = msk_roms.shape
etamax, ximax = sh2

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
 
zc=np.array([0])

zc=zc[::-1]

uavg=avgfile['u_eastward'][:]

vavg=avgfile['v_northward'][:]

time=avgfile['ocean_time'][:]/(24*60*60)

tempavg=avgfile['temp'][:]

uavg2=avgfile2['u_eastward'][:]

vavg2=avgfile2['v_northward'][:]

time2=avgfile2['ocean_time'][:]/(24*60*60)

tempavg2=avgfile2['temp'][:]

uavg=np.concatenate([uavg,uavg2], axis=0)
vavg=np.concatenate([vavg,vavg2], axis=0)
tempavg=np.concatenate([tempavg,tempavg2], axis=0)


time=np.concatenate([time,time2])

intu=np.zeros([uavg.shape[0],len(zc),x_roms.shape[0], x_roms.shape[1]])

intv=np.zeros([uavg.shape[0],len(zc),x_roms.shape[0], x_roms.shape[1]])

itemp=np.zeros([uavg.shape[0],len(zc),x_roms.shape[0], x_roms.shape[1]])

#zeta=avgfile['zeta'][:]

for j in range(intu.shape[2]):
  for k in range(intu.shape[3]):
    if (zr[-1,j,k] > zc.min()):
      zr[-1,j,k] = zc.min()

UNDEF=np.nan

for i in range(intu.shape[0]):
  for j in range(intu.shape[2]):
    for k in range(intu.shape[3]):
#      intu[i,:,j,k] = np.interp(-zc, -zr[:,j,k], uavg[i,:,j,k], right=UNDEF, left=UNDEF)
#      intv[i,:,j,k] = np.interp(-zc, -zr[:,j,k], vavg[i,:,j,k], right=UNDEF, left=UNDEF)
      itemp[i,:,j,k] = np.interp(-zc, -zr[:,j,k], tempavg[i,:,j,k], right=UNDEF, left=UNDEF)

intu=np.squeeze(intu)
intv=np.squeeze(intv)
itemp=np.squeeze(itemp)


begindate=avgfile['ocean_time'].units[14:]
begindate=dates.datestr2num(begindate)
romstime=begindate+time

########ostia

ostiafile=Dataset('METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2_1590691999912.nc')

ostialat=ostiafile['lat'][:]
ostialon=ostiafile['lon'][:]
ostiatemp=ostiafile['analysed_sst'][:]
ostiatemp=ostiatemp-273.15

timeostia=ostiafile['time'][:]

begindate=ostiafile['time'].units[14:]
begindate=dates.datestr2num(begindate)

timeostia=begindate + (timeostia/(24*60*60))

############################################################

intunewtime=np.zeros([len(timeostia),x_roms.shape[0], x_roms.shape[1]])

intvnewtime=np.zeros([len(timeostia),x_roms.shape[0], x_roms.shape[1]])

tempnewtime=np.zeros([len(timeostia),x_roms.shape[0], x_roms.shape[1]])


for j in range(intu.shape[1]):
  for k in range(intu.shape[2]):
#    intunewtime[:,j,k] = np.interp(timeostia, romstime, intu[:,j,k], right=UNDEF, left=UNDEF)
#    intvnewtime[:,j,k] = np.interp(timeostia, romstime, intv[:,j,k], right=UNDEF, left=UNDEF)
    tempnewtime[:,j,k] = np.interp(timeostia, romstime, itemp[:,j,k], right=UNDEF, left=UNDEF)
 
 
 
tempinterpostia=np.zeros(tempnewtime.shape)

from scipy.interpolate import interp2d

for i in range(tempinterpostia.shape[0]):
  tempinput = interp2d(ostialon, ostialat,ostiatemp[i,:])  
  tempinterpostia[i,:] = tempinput(x_roms[0,:], y_roms[:,0])
  

tempnewtime=np.ma.masked_where(tempnewtime>100, tempnewtime)
tempinterpostia=np.ma.masked_where(tempinterpostia>100, tempinterpostia)
tempinterpostia=np.ma.masked_where(tempinterpostia<-100, tempinterpostia)


from sklearn.metrics import mean_squared_error

temprmse=np.zeros(x_roms.shape)

for j in range(intu.shape[1]):
  for k in range(intu.shape[2]):
    temprmse[j,k] = np.sqrt(mean_squared_error(tempinterpostia[:,j,k],tempnewtime[:,j,k]))

temprmse = np.ma.masked_where(temprmse>100, temprmse)

a.min()

plt.pcolor(b[20,:],vmin=19,vmax=29);plt.colorbar();plt.show()
plt.pcolor(ostiatemp[20,:],vmin=19,vmax=29);plt.colorbar();plt.show()

plt.pcolor(temprmse);plt.colorbar();plt.show()


###############MUR
#source: https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41.html

murfile=Dataset('jplMURSST41_b166_d721_4325.nc')

murlat=murfile['latitude'][:]
murlon=murfile['longitude'][:]
murtemp=murfile['analysed_sst'][:]

timemur=murfile['time'][:]

begindate=murfile['time'].units[14:]
begindate=dates.datestr2num(begindate)

timemur=begindate + (timemur/(24*60*60))

intunewtime=np.zeros([len(timemur),x_roms.shape[0], x_roms.shape[1]])

intvnewtime=np.zeros([len(timemur),x_roms.shape[0], x_roms.shape[1]])

tempnewtime=np.zeros([len(timemur),x_roms.shape[0], x_roms.shape[1]])


for j in range(intu.shape[1]):
  for k in range(intu.shape[2]):
#    intunewtime[:,j,k] = np.interp(timemur, romstime, intu[:,j,k], right=UNDEF, left=UNDEF)
#    intvnewtime[:,j,k] = np.interp(timemur, romstime, intv[:,j,k], right=UNDEF, left=UNDEF)
    tempnewtime[:,j,k] = np.interp(timemur, romstime, itemp[:,j,k], right=UNDEF, left=UNDEF)
    

tempinterpmur=np.zeros(tempnewtime.shape)

from scipy.interpolate import interp2d

murtemp = np.ma.filled(murtemp,fill_value=100000000)
for i in range(tempinterpmur.shape[0]):
  tempinput = interp2d(murlon, murlat,murtemp[i,:])  
  tempinterpmur[i,:] = tempinput(x_roms[0,:], y_roms[:,0])

tempnewtime=np.ma.masked_where(tempnewtime>100, tempnewtime)
tempinterpmur=np.ma.masked_where(tempinterpmur>100, tempinterpmur)

#tempinterpmur=np.ma.masked_where(np.repeat(msk_roms[np.newaxis,:],len(timemur),axis=0)==0, tempinterpmur)

tempnewtimeNa = tempnewtime.copy()

#tempdif=tempnewtime-tempinterpmur 

tempdifNA=tempnewtimeNa-tempinterpmur 


plt.pcolor(np.squeeze(tempnewtime - tempinterpmur), vmax=2, vmin=-2);plt.colorbar();plt.show()  
  
  
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error

temprmsemur=np.zeros(x_roms.shape)
tempmaemur=np.zeros(x_roms.shape)


tempdif=tempnewtime-tempinterpmur 


for j in range(intu.shape[1]):
  for k in range(intu.shape[2]):
    temprmsemur[j,k] = np.sqrt(mean_squared_error(tempinterpmur[:,j,k],tempnewtime[:,j,k]))
    tempmaemur[j,k] = mean_absolute_error(tempinterpmur[:,j,k],tempnewtime[:,j,k])


temprmsemur = np.ma.masked_where(temprmsemur>1.5, temprmsemur)
tempmaemur = np.ma.masked_where(tempmaemur>1.5, tempmaemur)


plt.pcolor(temprmsemur,vmax=1.5);plt.colorbar();plt.show()

plt.pcolor(tempmaemur);plt.colorbar();plt.show()


plt.pcolor(tempdif[-1,:], vmin=-3, vmax=3);plt.colorbar();plt.show()

plt.pcolor(temprmse, vmax=1.5);plt.colorbar();plt.show()

peres, 2017

POES 1 day supper obs


from mpl_toolkits.basemap import Basemap

nlim=y_roms.max()
slim=y_roms.min()
wlim=x_roms.min()
elim=x_roms.max()

fig, axs = plt.subplots(figsize=(8,8))
map = Basemap(projection='merc', llcrnrlat=slim, urcrnrlat=nlim,llcrnrlon=wlim, urcrnrlon=elim,resolution='l')
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.fillcontinents(color='#ddaa66',lake_color='aqua')
parallels = np.arange(-30,-10,1)
map.drawparallels(parallels,labels=[1,0,0,1], linewidth=0.0)
meridians = np.arange(-55,-30,1.5)
map.drawmeridians(meridians,labels=[1,0,0,1], linewidth=0.0)
map.readshapefile('/mnt/c/Users/Fernando/Desktop/shape_unid_fed/lim_unidade_federacao_a', 'lim_unidade_federacao_a')
x, y = map(x_roms,y_roms)
#a=axs.pcolor(x,y, np.squeeze(temprmsemur),vmax=1.5)
a=axs.pcolor(x,y, np.squeeze(temprmsemur),vmax=1.5, cmap=plt.get_cmap('GnBu')) #source https://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3   sequential
divider = make_axes_locatable(axs)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar=fig.colorbar(a, shrink=0.8, extend='both', ax=axs, cax=cax)
cbar.ax.set_ylabel('RMSE Temperature ($^\circ$C)', rotation=270)
cbar.ax.get_yaxis().labelpad = 20
plt.show()



nlim=y_roms.max()
slim=y_roms.min()
wlim=x_roms.min()
elim=x_roms.max()

fig, axs = plt.subplots(figsize=(8,8))
map = Basemap(projection='merc', llcrnrlat=slim, urcrnrlat=nlim,llcrnrlon=wlim, urcrnrlon=elim,resolution='l')
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.fillcontinents(color='#ddaa66',lake_color='aqua')
parallels = np.arange(-30,-10,1)
map.drawparallels(parallels,labels=[1,0,0,1], linewidth=0.0)
meridians = np.arange(-55,-30,1.5)
map.drawmeridians(meridians,labels=[1,0,0,1], linewidth=0.0)
map.readshapefile('/mnt/c/Users/Fernando/Desktop/shape_unid_fed/lim_unidade_federacao_a', 'lim_unidade_federacao_a')
x, y = map(x_roms,y_roms)
#a=axs.pcolor(x,y, np.squeeze(temprmsemur),vmax=1.5)
a=axs.pcolor(x,y, np.squeeze(tempdif),vmax=1.5, vmin=-1.5, cmap=plt.get_cmap('seismic')) #source https://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3   sequential
divider = make_axes_locatable(axs)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar=fig.colorbar(a, shrink=0.8, extend='both', ax=axs, cax=cax)
cbar.ax.set_ylabel('RMSE Temperature ($^\circ$C)', rotation=270)
cbar.ax.get_yaxis().labelpad = 20
plt.show()


###SUBPLOT 2 GRIDS
from mpl_toolkits.basemap import Basemap


nlim=y_roms.max()
slim=y_roms.min()
wlim=x_roms.min()
elim=x_roms.max()

mint=-2
maxt=2

fig, ax1 = plt.subplots(2,1, figsize=(8,15))
map = Basemap(projection='merc', llcrnrlat=slim, urcrnrlat=nlim,llcrnrlon=wlim, urcrnrlon=elim,resolution='l', ax=ax1[0])
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.fillcontinents(color='#ddaa66',lake_color='aqua')
parallels = np.arange(-30,-10,2)
map.drawparallels(parallels,labels=[1,0,0,1], linewidth=0.0, fontsize=9)
meridians = np.arange(-55,-20,3)
map.drawmeridians(meridians,labels=[1,0,0,1], linewidth=0.0, fontsize=9)
map.readshapefile('/mnt/c/Users/Fernando/Desktop/shape_unid_fed/lim_unidade_federacao_a', 'lim_unidade_federacao_a')
x, y = map(x_roms,y_roms)
aa=ax1[0].pcolor(x,y, np.squeeze(tempdif),vmax=maxt, vmin=mint, cmap=plt.get_cmap('seismic')) #source https://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3   sequential
ax1[0].title.set_text('D.A. RUN')
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=12)



map = Basemap(projection='merc', llcrnrlat=slim, urcrnrlat=nlim,llcrnrlon=wlim, urcrnrlon=elim,resolution='l', ax=ax1[1])
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.fillcontinents(color='#ddaa66',lake_color='aqua')
parallels = np.arange(-30,-10,2)
map.drawparallels(parallels,labels=[1,0,0,1], linewidth=0.0, fontsize=9)
meridians = np.arange(-55,-20,3)
map.drawmeridians(meridians,labels=[1,0,0,1], linewidth=0.0, fontsize=9)
map.readshapefile('/mnt/c/Users/Fernando/Desktop/shape_unid_fed/lim_unidade_federacao_a', 'lim_unidade_federacao_a')
x, y = map(x_roms,y_roms)
a=ax1[1].pcolor(x,y, np.squeeze(tempdifNA),vmax=maxt, vmin=mint, cmap=plt.get_cmap('seismic')) #source https://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3   sequential
ax1[1].title.set_text('FREE RUN')
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=12)




p0 = ax1[0].get_position().get_points().flatten()
p1 = ax1[1].get_position().get_points().flatten()

ax_cbar = fig.add_axes([0.78, 0.3, 0.03, 0.4])

cbar = plt.colorbar(aa, cax=ax_cbar, extend='both')

cbar.ax.set_ylabel('Temperature ($^\circ$C)', rotation=270)
cbar.ax.get_yaxis().labelpad = 11
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=12)
text.set_font_properties(font)

plt.show()
