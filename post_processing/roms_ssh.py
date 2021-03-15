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
avgfile=Dataset('SSH_Media_BRSE_fwd_outer1.nc')
#avgfile=Dataset('his_no_ass.nc')
#avgfile2=Dataset('HIS_FILE_20200421_5D0-20200428_5D0_hind_correct_year_WEAK_menor_azul_nopline_0007.nc')

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

zetaavg=avgfile['zeta'][:]

time=avgfile['ocean_time'][:]/(24*60*60)

tempavg=avgfile['temp'][:]

uavg2=avgfile2['u_eastward'][:]

vavg2=avgfile2['v_northward'][:]

time2=avgfile2['ocean_time'][:]/(24*60*60)

tempavg2=avgfile2['temp'][:]

zetaavg2=avgfile2['zeta'][:]


uavg=np.concatenate([uavg,uavg2], axis=0)
vavg=np.concatenate([vavg,vavg2], axis=0)
tempavg=np.concatenate([tempavg,tempavg2], axis=0)

zetaavg=np.concatenate([zetaavg,zetaavg2], axis=0)

time=np.concatenate([time[:],time2[:]])

izeta=np.squeeze(zetaavg)


begindate=avgfile['ocean_time'].units[14:]
begindate=dates.datestr2num(begindate)
romstime=begindate+time


izeta=np.ma.masked_where(izeta>100, izeta)

tim=-1
a=plt.pcolor(x_roms,y_roms, np.squeeze(izeta[tim,:]),  vmin=-0.1, vmax=0.4)
plt.text(-46,-19,dates.num2date(romstime[tim]).strftime("%Y%m%d - %-I %p"), fontsize=12, fontweight='bold',
        bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 5})
plt.colorbar(a)
plt.show()


########ostia

#ostiafile=Dataset('dataset-duacs-rep-global-merged-allsat-phy-l4_1591217004924.nc')
#ostiafile=Dataset('dataset-duacs-nrt-global-merged-allsat-phy-l4_1605719127190.nc')
ostiafile=Dataset('alt_1.nc')

ostialat=ostiafile['latitude'][:]
#ostialon=ostiafile['longitude'][:] - 360
ostialon=ostiafile['longitude'][:]
ostiatemp=ostiafile['adt'][:]

#anoostiatemp=ostiafile['sla'][:]

#ostiatemp = ostiatemp - anoostiatemp
ostiatemp = ostiatemp - 0.6719 + 0.2187

timeostia=ostiafile['time'][:]

begindate=ostiafile['time'].units[11:]
begindate=dates.datestr2num(begindate)

timeostia=begindate + (timeostia)


############################################################

intunewtime=np.zeros([len(timeostia),x_roms.shape[0], x_roms.shape[1]])

intvnewtime=np.zeros([len(timeostia),x_roms.shape[0], x_roms.shape[1]])

tempnewtime=np.zeros([len(timeostia),x_roms.shape[0], x_roms.shape[1]])

UNDEF=np.nan

for j in range(izeta.shape[1]):
  for k in range(izeta.shape[2]):
    tempnewtime[:,j,k] = np.interp(timeostia, romstime, izeta[:,j,k], right=UNDEF, left=UNDEF)

from scipy.interpolate import interp2d


tempinterpostia=np.zeros(tempnewtime.shape)


for i in range(tempinterpostia.shape[0]):
  tempinput = interp2d(ostialon, ostialat,ostiatemp[i,:])  
  tempinterpostia[i,:] = tempinput(x_roms[0,:], y_roms[:,0])


tempinterpostia=np.ma.masked_where(tempinterpostia<-100, tempinterpostia)
tempnewtime=np.ma.masked_where(tempnewtime>100, tempnewtime)


from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error

temprmsemur=np.zeros(x_roms.shape)
tempmaemur=np.zeros(x_roms.shape)



for j in range(izeta.shape[1]):
  for k in range(izeta.shape[2]):
    temprmsemur[j,k] = np.sqrt(mean_squared_error(tempinterpostia[:,j,k],tempnewtime[:,j,k]))
    tempmaemur[j,k] = mean_absolute_error(tempinterpostia[:,j,k],tempnewtime[:,j,k])
   

   
temprmsemur = np.ma.masked_where(temprmsemur>3., temprmsemur)
tempmaemur = np.ma.masked_where(tempmaemur>3., tempmaemur)



tempnewtimeNa = tempnewtime.copy()



tempdif=tempnewtime-tempinterpostia 

tempdifNA=tempnewtimeNa-tempinterpostia 


plt.pcolor(x_roms, y_roms, np.squeeze(tempinterpostia));plt.colorbar();plt.show()


plt.pcolor(np.squeeze(tempnewtime));plt.colorbar();plt.show()


plt.pcolor(np.squeeze(tempnewtimeNa));plt.colorbar();plt.show()



plt.pcolor(np.squeeze(tempinterpostia));plt.colorbar();plt.show()


plt.pcolor(np.squeeze(tempnewtime)+0.5);plt.colorbar();plt.show()

 

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
a=axs.pcolor(x,y, np.squeeze(tempmaemur))
#a=axs.pcolor(x,y, np.squeeze(tempdif),cmap=plt.get_cmap('seismic')) #source https://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3   sequential
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

mint=-0.1
maxt=0.1

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
#aa=ax1[0].pcolor(x,y, np.squeeze(tempdif),cmap=plt.get_cmap('seismic')) #source https://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3   sequential
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
#a=ax1[1].pcolor(x,y, np.squeeze(tempdifNA),cmap=plt.get_cmap('seismic')) #source https://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3   sequential
ax1[1].title.set_text('FREE RUN')
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=12)




p0 = ax1[0].get_position().get_points().flatten()
p1 = ax1[1].get_position().get_points().flatten()

ax_cbar = fig.add_axes([0.78, 0.3, 0.03, 0.4])

cbar = plt.colorbar(aa, cax=ax_cbar, extend='both')

cbar.ax.set_ylabel('Altimetry (m)', rotation=270)
cbar.ax.get_yaxis().labelpad = 11
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=12)
text.set_font_properties(font)

plt.show()