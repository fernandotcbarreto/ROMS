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

lista = sorted(glob.glob('ocean_BRSE_his_b5_000[0-9]*'))

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

uavg=np.array(avgfile.u.isel(s_rho=lst).values)


uavg = 0.5*(uavg[:,:,:,1:]+uavg[:,:,:,:-1])
lonu = 0.5*(lonu[:,1:]+lonu[:,:-1])
latu = 0.5*(latu[:,1:]+latu[:,:-1])


uavg=uavg[:,:,1:-1,:]
lonu=lonu[1:-1,:]
latu=latu[1:-1,:]

##np.array(avgfile.u.isel(s_rho=lst).values)

vavg=np.array(avgfile.v.isel(s_rho=lst).values)

vavg = 0.5*(vavg[:,:,1:,:]+vavg[:,:,:-1,:])
latv = 0.5*(latv[1:,:]+latv[:-1,:])
lonv = 0.5*(lonv[1:,:]+lonv[:-1,:])

vavg=vavg[:,:,:,1:-1]
lonv=lonv[:,1:-1]
latv=latv[:,1:-1]


tempavg=np.array(avgfile.temp.isel(s_rho=lst).values)
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


lon  = x_roms
lat  = y_roms
u    = intu[-1,:]
v    = intv[-1,:]
t    = itemp[-1,:]

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




for l in range(-1):
  count += 1
  bc =  ( (g*alpha)/theta_z ) * ( unot[l]*tnot[l]*dTdx + vnot[l]*tnot[l]*dTdy )

plt.pcolor(x_roms,y_roms,bc);plt.colorbar();plt.show()

