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

avgfile=Dataset('HIS_FILE_20200421_5D0-20200428_5D0_hind_correct_year_WEAK_menor_azul_nopline_0009.nc')
avgfile2=Dataset('HIS_FILE_20200421_5D0-20200428_5D0_hind_correct_year_WEAK_menor_azul_nopline_0010.nc')

fname_grd = 'azul_grd_era_NEW_menor_azul.nc'  ## Smooth topography grid.   


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

time=avgfile['ocean_time'][:]

tempavg=avgfile['temp'][:]

uavg2=avgfile2['u_eastward'][:]

vavg2=avgfile2['v_northward'][:]

time2=avgfile2['ocean_time'][:]

tempavg2=avgfile2['temp'][:]

zetaavg2=avgfile2['zeta'][:]


uavg=np.concatenate([uavg,uavg2], axis=0)
vavg=np.concatenate([vavg,vavg2], axis=0)
tempavg=np.concatenate([tempavg,tempavg2], axis=0)

zetaavg=np.concatenate([zetaavg,zetaavg2], axis=0)

time=np.concatenate([time[:]/(24*60*60),time2[:]/(24*60*60)])

izeta=np.squeeze(zetaavg)


begindate=avgfile['ocean_time'].units[14:]
begindate=dates.datestr2num(begindate)
romstime=begindate+time


izeta=np.ma.masked_where(izeta>100, izeta)

tim=86
a=plt.pcolor(x_roms,y_roms, np.squeeze(izeta[tim,:]),  vmin=-0.1, vmax=0.4)
plt.text(-46,-19,dates.num2date(romstime[tim]).strftime("%m%d - %-I %p"), fontsize=12, fontweight='bold',
        bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 5})
plt.colorbar(a)
plt.show()


########ostia

ostiafile=Dataset('dataset-duacs-rep-global-merged-allsat-phy-l4_1591217004924.nc')

ostialat=ostiafile['latitude'][:]
ostialon=ostiafile['longitude'][:] - 360
ostiatemp=ostiafile['adt'][:]



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


plt.pcolor(x_roms, y_roms, np.squeeze(tempinterpostia));plt.colorbar();plt.show()


plt.pcolor(np.squeeze(tempnewtime));plt.colorbar();plt.show()





plt.pcolor(np.squeeze(tempinterpostia));plt.colorbar();plt.show()


plt.pcolor(np.squeeze(tempnewtime)+0.5);plt.colorbar();plt.show()

 
