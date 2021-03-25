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
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

#lista = sorted(glob.glob('ocean_BRSE_his_b5_000[0-9]*'))

lista = sorted(glob.glob('/mnt/share/Modelos/BRSE_2014_2016/RESULTADOS/ocean_BRSE_his_b1_00[5-9][0-9]*'))

avgfile = xr.open_mfdataset(lista, concat_dim='ocean_time')

hh=np.array(avgfile['ocean_time'][:].dt.strftime('%Y-%m-%d  %H:00:00'))

hhdates=dates.datestr2num(hh)

fname_grd = 'BRSE_2012_GRD.nc'     #CF 1/36


## Load ROMS grid.
grd = Dataset(fname_grd)
x_roms = grd.variables['lon_rho'][:]
y_roms = grd.variables['lat_rho'][:]
msk_roms = grd.variables['mask_rho'][:]
msk_romsv = grd.variables['mask_v'][:]
h_roms = grd.variables['h'][:]
sh2 = msk_roms.shape
etamax, ximax = sh2

x_roms=x_roms[1:-1,1:-1]
y_roms=y_roms[1:-1,1:-1]
####h_roms=h_roms[1:-1,1:-1]



theta_b = 0.4
theta_s = 5.0
tcline = 3.
klevels = 30
Vtransform = 1
Vstretching = 1
Spherical = True

lst=list(np.arange(0,30))

if Vstretching==4:
  scoord = s_coordinate_4(h_roms, theta_b, theta_s, tcline, klevels)   #zeta is not used in the computation 
elif Vstretching==2:
  scoord = s_coordinate_2(h_roms, theta_b, theta_s, tcline, klevels)
elif Vstretching==1:
  scoord = s_coordinate(h_roms, theta_b, theta_s, tcline, klevels)

zr = -scoord.z_r[:]

zr=zr[lst,1:-1,1:-1] 

zc=np.array([4.940250e-01, 1.541375e+00, 2.645669e+00, 3.819495e+00,                                                                 5.078224e+00, 6.440614e+00, 7.929560e+00, 9.572997e+00,                                                                 1.140500e+01, 1.346714e+01, 1.581007e+01, 1.849556e+01,                                                                 2.159882e+01, 2.521141e+01, 2.944473e+01, 3.443415e+01,                                                                 4.034405e+01, 4.737369e+01, 5.576429e+01, 6.580727e+01,                                                                 7.785385e+01, 9.232607e+01, 1.097293e+02, 1.306660e+02,                                                                 1.558507e+02, 1.861256e+02, 2.224752e+02, 2.660403e+02,                                                                 3.181274e+02, 3.802130e+02, 4.539377e+02, 5.410889e+02,                                                                 6.435668e+02, 7.633331e+02, 9.023393e+02, 1.062440e+03,                                                                 1.245291e+03, 1.452251e+03, 1.684284e+03, 1.941893e+03,                                                                 2.225078e+03, 2.533336e+03, 2.865703e+03, 3.220820e+03,                                                                 3.597032e+03, 3.992484e+03, 4.405224e+03, 4.833291e+03,                                                               5.274784e+03, 5.727917e+03])



zc=zc[::-1]


latv=np.array(avgfile['lat_v'][:])
lonv=np.array(avgfile['lon_v'][:])
lonu=np.array(avgfile['lon_u'][:])
latu=np.array(avgfile['lat_u'][:])

#uavg=np.array(avgfile.u.isel(s_rho=lst).values)
uavg=np.array(avgfile.u.isel(s_rho=lst).groupby('ocean_time.month').mean().values)


uavg = 0.5*(uavg[:,:,:,1:]+uavg[:,:,:,:-1])
lonu = 0.5*(lonu[:,1:]+lonu[:,:-1])
latu = 0.5*(latu[:,1:]+latu[:,:-1])


uavg=uavg[:,:,1:-1,:]
lonu=lonu[1:-1,:]
latu=latu[1:-1,:]

##np.array(avgfile.u.isel(s_rho=lst).values)

#vavg=np.array(avgfile.v.isel(s_rho=lst).values)
vavg=np.array(avgfile.v.isel(s_rho=lst).groupby('ocean_time.month').mean().values)


vavg = 0.5*(vavg[:,:,1:,:]+vavg[:,:,:-1,:])
latv = 0.5*(latv[1:,:]+latv[:-1,:])
lonv = 0.5*(lonv[1:,:]+lonv[:-1,:])

vavg=vavg[:,:,:,1:-1]
lonv=lonv[:,1:-1]
latv=latv[:,1:-1]


#tempavg=np.array(avgfile.temp.isel(s_rho=lst).values)

tempavg=np.array(avgfile.temp.isel(s_rho=lst).groupby('ocean_time.month').mean().values)


tempavg=tempavg[:,:,1:-1,1:-1]


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
      intu[i,:,j,k] = np.interp(-zc, -zr[:,j,k], uavg[i,:,j,k], right=UNDEF, left=UNDEF)
      intv[i,:,j,k] = np.interp(-zc, -zr[:,j,k], vavg[i,:,j,k], right=UNDEF, left=UNDEF)
      itemp[i,:,j,k] = np.interp(-zc, -zr[:,j,k], tempavg[i,:,j,k], right=UNDEF, left=UNDEF)




import numpy as np
import sys
import scipy.io as sp
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from mpl_toolkits.basemap import Basemap
import datetime as dt
import seawater as sw
import netCDF4 as nc
from matplotlib.patches import Polygon
from cookb_signalsmooth import smooth

tim=0
lon  = x_roms
lat  = y_roms
u    = intu[tim,:]
v    = intv[tim,:]
t    = itemp[tim,:]

zlev=0

u = np.ma.masked_where(np.abs(u) > 2, u)
v = np.ma.masked_where(np.abs(v) > 2, v)
t = np.ma.masked_where(np.abs(v) > 30, t)

u = np.ma.masked_invalid(u)
v = np.ma.masked_invalid(v)
t = np.ma.masked_invalid(t)


lm, jm, im = u.shape

# computing averages and perturbations ####################################
ubar, vbar, tbar = u.mean(axis=0), v.mean(axis=0), t.mean(axis=0)
unot, vnot, tnot = u*0, v*0, t*0

for l in range(lm):
  unot[l,...] = u[l,...] - ubar
  vnot[l,...] = v[l,...] - vbar
  tnot[l,...] = t[l,...] - tbar

# baroclinic energy  conversions ##########################################
bec = []
g = 9.8

theta_z = tbar.mean()
alpha = sw.alpha(35, theta_z, zlev)
gradT = np.gradient(tbar)
dTdx  = gradT[1] / (0.05 * 111000)
dTdy  = gradT[0] / (0.05 * 111000)

count = 0

for l in range(lm):
  count += 1
  bc =  ( (g*alpha)/theta_z ) * ( unot[l]*tnot[l]*dTdx + vnot[l]*tnot[l]*dTdy )
  bec.append( bc.mean() )

bec = np.array(bec)

bec = smooth(bec, window_len=10, window='hanning')



sp=4
#ly=16
ly=-1
for l in range(ly):
  count += 1
  bc =  ( (g*alpha)/theta_z ) * ( unot[l]*tnot[l]*dTdx + vnot[l]*tnot[l]*dTdy )

#plt.pcolor(x_roms,y_roms,bc*10**(7), cmap=plt.get_cmap('seismic'), vmin=-0.5, vmax=0.5);plt.colorbar();
#plt.quiver(x_roms[0:-1:sp, 0:-1:sp],y_roms[0:-1:sp, 0:-1:sp],u[ly,0:-1:sp, 0:-1:sp], v[ly,0:-1:sp, 0:-1:sp])
#plt.show()


lon1=x_roms.copy()
lat1=y_roms.copy()

for tim in [tim]:
  sp=6
  up=u[ly,:,:]
  vp=v[ly,:,:]
  uu=up
  vv=vp
#  val=np.sqrt((uu**2)+(vv**2))
#  begindate=avgfile['ocean_time'].units[14:]
#  begindate=dates.datestr2num(begindate)
#  figdates=dates.num2date(begindate+time)
  nlim=lat1.max()
  slim=lat1.min()
  wlim=lon1.min()
  elim=lon1.max()
  fig, ax1 = plt.subplots(figsize=(8,8))
  map = Basemap(projection='cyl', llcrnrlat=slim, urcrnrlat=nlim,llcrnrlon=wlim, urcrnrlon=elim,resolution='l')
  map.drawcoastlines(linewidth=0.25)
  map.drawcountries(linewidth=0.25)
  map.fillcontinents(color='silver',lake_color='white')
  parallels = np.arange(-50,-9,3)
  map.drawparallels(parallels,labels=[1,0,0,1], linewidth=0.0)
  meridians = np.arange(-65,-15,4)
  map.drawmeridians(meridians,labels=[1,0,0,1], linewidth=0.0)
  map.readshapefile('/mnt/c/Users/Fernando/Desktop/shape_unid_fed/lim_unidade_federacao_a', 'lim_unidade_federacao_a')
#  map.drawmapscale(-39., -24.5, -39, -24.5, 100, barstyle='fancy', fontsize = 8, yoffset=5000)
  x_roms, y_roms = map(lon1,lat1)
#
  a=ax1.pcolor(x_roms,y_roms,bc*10**(7), cmap=plt.get_cmap('seismic'), vmin=-0.5, vmax=0.5)
  ax1.quiver(x_roms[0:-1:sp, 0:-1:sp],y_roms[0:-1:sp, 0:-1:sp],uu[0:-1:sp, 0:-1:sp], vv[0:-1:sp, 0:-1:sp], alpha=1)
#  ax1.contourf(x_roms,y_roms,h_roms,levels=[0,600],colors=('gainsboro'))
#  plt.text(-49.5,-17,figdates[tim].strftime("%Y/%b/%d - %-I %p"), fontsize=10, fontweight='bold',
#        bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 5})
  divider = make_axes_locatable(ax1)
  #cax = divider.append_axes("right", size="5%", pad=0.05)
  cax = fig.add_axes([0.91, 0.3, 0.02, 0.38])
  cbar=fig.colorbar(a, shrink=0.8, extend='both', ax=ax1, cax=cax)
  cbar.ax.set_ylabel('$m^2s^{-3}$ ($10^{-7}$)', rotation=270)
  cbar.ax.get_yaxis().labelpad = 11
  cbar.ax.tick_params(labelsize=9)
  text = cbar.ax.yaxis.label
  font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=11)
  text.set_font_properties(font)
#
#
#
#
  left, bottom, width, height = [0.55, 0.17, 0.35, 0.35]
  axins = fig.add_axes([left, bottom, width, height])
#
#
  nlim=-16
  slim=-24.5
  wlim=-43
  elim=-31
#
  map = Basemap(projection='cyl', llcrnrlat=slim, urcrnrlat=nlim,llcrnrlon=wlim, urcrnrlon=elim,resolution='l')
  map.drawcoastlines(linewidth=0.25)
  map.drawcountries(linewidth=0.25)
  map.fillcontinents(color='silver',lake_color='white')
  parallels = np.arange(-40,-15,1)
  map.drawparallels(parallels,labels=[0,0,0,0], linewidth=0.0)
  meridians = np.arange(-55,-20,1)
  map.drawmeridians(meridians,labels=[0,0,0,0], linewidth=0.0)
  map.readshapefile('/mnt/c/Users/Fernando/Desktop/shape_unid_fed/lim_unidade_federacao_a', 'lim_unidade_federacao_a')
#
#
#  axins.pcolor(x_roms2,y_roms2, np.squeeze(val),vmin=val.min(),vmax=1.0)
#  axins.contourf(x_roms2,y_roms2,h_roms2,levels=[0,20],colors=('gainsboro'))
#  axins.quiver(x_roms2[0:-1:sp, 0:-1:sp],y_roms2[0:-1:sp, 0:-1:sp],uu[0:-1:sp, 0:-1:sp], vv[0:-1:sp, 0:-1:sp])
  a=axins.pcolor(x_roms,y_roms,bc*10**(7), cmap=plt.get_cmap('seismic'), vmin=-0.5, vmax=0.5)
  sp=10
#  axins.quiver(x_roms[0:-1:sp, 0:-1:sp],y_roms[0:-1:sp, 0:-1:sp],uu[0:-1:sp, 0:-1:sp], vv[0:-1:sp, 0:-1:sp], alpha=1)
  axins.contourf(x_roms,y_roms,h_roms,levels=[0,600],colors=('gainsboro'))
#
  mark_inset(ax1, axins, loc1=1, loc2=3, fc="none", ec="black")
#
#

plt.savefig('EC_00.png', dpi=200, bbox_inches='tight', transparent=False)

plt.show()

