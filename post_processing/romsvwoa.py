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

########WOA TEMP


file=Dataset('woa18_decav_t00_04.nc')
lat=file['lat'][3].filled()
lon=file['lon'][3].filled()
tempwoa=file['t_an'][0,:,3,3]

lista = sorted(glob.glob('/mnt/share/Modelos/BRSE_2014_2016/RESULTADOS/ocean_BRSE_his_b1_0[0-9][0-9][0-9]*'))
#lista = sorted(glob.glob('/mnt/share/Modelos/BRSE_2014_2016/RESULTADOS/ocean_BRSE_his_b1_0001*'))


def select_point(ds):
    return ds.drop(['u','v','ubar','vbar','zeta'])


avgfile = xr.open_mfdataset(lista, concat_dim='ocean_time', preprocess=select_point)


hh=np.array(avgfile['ocean_time'][:].dt.strftime('%Y-%m-%d  %H:00:00'))

hhdates=dates.datestr2num(hh)

fname_grd = 'BRSE_2012_GRD.nc'     #CF 1/36

tlst=list(np.arange(0,hh.shape[0],1))
     

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

ktl=3

lglst=np.arange(lgmin-ktl, lgmin+ktl)
ltlst=np.arange(ltmin-ktl, ltmin+ktl)


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

zc=np.array([4.940250e-01, 1.541375e+00, 2.645669e+00, 3.819495e+00,                                                                 5.078224e+00, 6.440614e+00, 7.929560e+00, 9.572997e+00,                                                                 1.140500e+01, 1.346714e+01, 1.581007e+01, 1.849556e+01,                                                                 2.159882e+01, 2.521141e+01, 2.944473e+01, 3.443415e+01,                                                                 4.034405e+01, 4.737369e+01, 5.576429e+01, 6.580727e+01,                                                                 7.785385e+01, 9.232607e+01, 1.097293e+02, 1.306660e+02,                                                                 1.558507e+02, 1.861256e+02, 2.224752e+02, 2.660403e+02,                                                                 3.181274e+02, 3.802130e+02, 4.539377e+02, 5.410889e+02,                                                                 6.435668e+02, 7.633331e+02, 9.023393e+02, 1.062440e+03,                                                                 1.245291e+03, 1.452251e+03, 1.684284e+03, 1.941893e+03,                                                                 2.225078e+03, 2.533336e+03, 2.865703e+03, 3.220820e+03,                                                                 3.597032e+03, 3.992484e+03, 4.405224e+03, 4.833291e+03,                                                               5.274784e+03, 5.727917e+03])

zc=np.array([0.000000000000000000e+00,5.000000000000000000e+00,1.000000000000000000e+01,1.500000000000000000e+01,2.000000000000000000e+01,2.500000000000000000e+01,3.000000000000000000e+01,3.500000000000000000e+01,4.000000000000000000e+01,4.500000000000000000e+01,5.000000000000000000e+01,5.500000000000000000e+01,6.000000000000000000e+01,6.500000000000000000e+01,7.000000000000000000e+01,7.500000000000000000e+01,8.000000000000000000e+01,8.500000000000000000e+01,9.000000000000000000e+01,9.500000000000000000e+01,1.000000000000000000e+02,1.250000000000000000e+02,1.500000000000000000e+02,1.750000000000000000e+02,2.000000000000000000e+02,2.250000000000000000e+02,2.500000000000000000e+02,2.750000000000000000e+02,3.000000000000000000e+02,3.250000000000000000e+02,3.500000000000000000e+02,3.750000000000000000e+02,4.000000000000000000e+02,4.250000000000000000e+02,4.500000000000000000e+02,4.750000000000000000e+02,5.000000000000000000e+02,5.500000000000000000e+02,6.000000000000000000e+02,6.500000000000000000e+02,7.000000000000000000e+02,7.500000000000000000e+02,8.000000000000000000e+02,8.500000000000000000e+02,9.000000000000000000e+02,9.500000000000000000e+02,1.000000000000000000e+03,1.050000000000000000e+03,1.100000000000000000e+03,1.150000000000000000e+03,1.200000000000000000e+03,1.250000000000000000e+03,1.300000000000000000e+03,1.350000000000000000e+03,1.400000000000000000e+03,1.450000000000000000e+03,1.500000000000000000e+03,1.550000000000000000e+03,1.600000000000000000e+03,1.650000000000000000e+03,1.700000000000000000e+03,1.750000000000000000e+03,1.800000000000000000e+03,1.850000000000000000e+03,1.900000000000000000e+03,1.950000000000000000e+03,2.000000000000000000e+03,2.100000000000000000e+03,2.200000000000000000e+03,2.300000000000000000e+03,2.400000000000000000e+03,2.500000000000000000e+03,2.600000000000000000e+03,2.700000000000000000e+03,2.800000000000000000e+03,2.900000000000000000e+03,3.000000000000000000e+03,3.100000000000000000e+03,3.200000000000000000e+03,3.300000000000000000e+03,3.400000000000000000e+03,3.500000000000000000e+03,3.600000000000000000e+03,3.700000000000000000e+03,3.800000000000000000e+03,3.900000000000000000e+03,4.000000000000000000e+03,4.100000000000000000e+03,4.200000000000000000e+03,4.300000000000000000e+03,4.400000000000000000e+03,4.500000000000000000e+03,4.600000000000000000e+03,4.700000000000000000e+03,4.800000000000000000e+03,4.900000000000000000e+03,5.000000000000000000e+03,5.100000000000000000e+03,5.200000000000000000e+03,5.300000000000000000e+03,5.400000000000000000e+03,5.500000000000000000e+03])
zc=zc[::-1]


tempavg=np.array(avgfile.temp.isel(s_rho=lst).isel(eta_rho=ltlst).isel(xi_rho=lglst).groupby('ocean_time.month').mean().values)

a,b=np.meshgrid(lglst, ltlst) 

tempavg=tempavg.mean(axis=0)


for j in range(x_roms.shape[0]):
  for k in range(x_roms.shape[1]):
    if (zr[-1,j,k] > zc.min()):
      zr[-1,j,k] = zc.min()

zr=zr[:,ltmin-ktl:ltmin+ktl, lgmin-ktl:lgmin+ktl]
x_roms=x_roms[ltmin-ktl:ltmin+ktl, lgmin-ktl:lgmin+ktl]
y_roms=y_roms[ltmin-ktl:ltmin+ktl, lgmin-ktl:lgmin+ktl]

itemp=np.zeros([len(zc),x_roms.shape[0], x_roms.shape[1]])


UNDEF=np.nan

for j in range(itemp.shape[1]):
  for k in range(itemp.shape[2]):
    itemp[:,j,k] = np.interp(-zc, -zr[:,j,k], tempavg[:,j,k], right=UNDEF, left=UNDEF)

itemp=np.squeeze(itemp)

temproms=np.zeros([len(zc)])

for i in range(len(vsitu)):
  temproms[i]=griddata((x_roms.ravel(),y_roms.ravel()), itemp[i,:].ravel(), (lon,lat))


temproms=temproms[::-1]



tempwoa=tempwoa.filled() 
rnan=np.where(isnan(temproms))
tempwoa[rnan]=np.nan
tempwoa[tempwoa>100]=np.nan


plt.close()
plt.style.use('ggplot')
fig, ax1 = plt.subplots(figsize=(5,7))
ax1.plot(temproms,-file['depth'][:], linewidth=2,color='blue', label='ROMS')
ax1.plot(tempwoa,-file['depth'][:], linewidth=2, color='red', label='WOA')
ax1.set_xlabel('Temperatura ($^\circ$C)', fontsize=14)
ax1.set_ylabel('Profundidade (m)',fontsize=14)
ax1.grid(alpha=0.5)
legend=ax1.legend(loc=4, fancybox=True, framealpha=1, shadow=True, borderpad=1, fontsize='small')
ax1.set_ylim([-4200,0])
plt.savefig('woa_temp.png', dpi=200, bbox_inches='tight', linestyle='.', transparent=False)




########WOA SALT



file=Dataset('woa18_decav_s00_04.nc')
lat=file['lat'][3].filled()
lon=file['lon'][3].filled()
saltwoa=file['s_an'][0,:,3,3]

#lista = sorted(glob.glob('/mnt/share/Modelos/BRSE_2014_2016/RESULTADOS/ocean_BRSE_his_b1_0[0-9][0-9][0-9]*'))
lista = sorted(glob.glob('/mnt/share/Modelos/BRSE_2014_2016/RESULTADOS/ocean_BRSE_his_b1_0001*'))


def select_point(ds):
    return ds.drop(['u','v','ubar','vbar','zeta'])


avgfile = xr.open_mfdataset(lista, concat_dim='ocean_time', preprocess=select_point)


hh=np.array(avgfile['ocean_time'][:].dt.strftime('%Y-%m-%d  %H:00:00'))

hhdates=dates.datestr2num(hh)

fname_grd = 'BRSE_2012_GRD.nc'     #CF 1/36

tlst=list(np.arange(0,hh.shape[0],1))
     

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

ktl=3

lglst=np.arange(lgmin-ktl, lgmin+ktl)
ltlst=np.arange(ltmin-ktl, ltmin+ktl)


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

zc=np.array([4.940250e-01, 1.541375e+00, 2.645669e+00, 3.819495e+00,                                                                 5.078224e+00, 6.440614e+00, 7.929560e+00, 9.572997e+00,                                                                 1.140500e+01, 1.346714e+01, 1.581007e+01, 1.849556e+01,                                                                 2.159882e+01, 2.521141e+01, 2.944473e+01, 3.443415e+01,                                                                 4.034405e+01, 4.737369e+01, 5.576429e+01, 6.580727e+01,                                                                 7.785385e+01, 9.232607e+01, 1.097293e+02, 1.306660e+02,                                                                 1.558507e+02, 1.861256e+02, 2.224752e+02, 2.660403e+02,                                                                 3.181274e+02, 3.802130e+02, 4.539377e+02, 5.410889e+02,                                                                 6.435668e+02, 7.633331e+02, 9.023393e+02, 1.062440e+03,                                                                 1.245291e+03, 1.452251e+03, 1.684284e+03, 1.941893e+03,                                                                 2.225078e+03, 2.533336e+03, 2.865703e+03, 3.220820e+03,                                                                 3.597032e+03, 3.992484e+03, 4.405224e+03, 4.833291e+03,                                                               5.274784e+03, 5.727917e+03])

zc=np.array([0.000000000000000000e+00,5.000000000000000000e+00,1.000000000000000000e+01,1.500000000000000000e+01,2.000000000000000000e+01,2.500000000000000000e+01,3.000000000000000000e+01,3.500000000000000000e+01,4.000000000000000000e+01,4.500000000000000000e+01,5.000000000000000000e+01,5.500000000000000000e+01,6.000000000000000000e+01,6.500000000000000000e+01,7.000000000000000000e+01,7.500000000000000000e+01,8.000000000000000000e+01,8.500000000000000000e+01,9.000000000000000000e+01,9.500000000000000000e+01,1.000000000000000000e+02,1.250000000000000000e+02,1.500000000000000000e+02,1.750000000000000000e+02,2.000000000000000000e+02,2.250000000000000000e+02,2.500000000000000000e+02,2.750000000000000000e+02,3.000000000000000000e+02,3.250000000000000000e+02,3.500000000000000000e+02,3.750000000000000000e+02,4.000000000000000000e+02,4.250000000000000000e+02,4.500000000000000000e+02,4.750000000000000000e+02,5.000000000000000000e+02,5.500000000000000000e+02,6.000000000000000000e+02,6.500000000000000000e+02,7.000000000000000000e+02,7.500000000000000000e+02,8.000000000000000000e+02,8.500000000000000000e+02,9.000000000000000000e+02,9.500000000000000000e+02,1.000000000000000000e+03,1.050000000000000000e+03,1.100000000000000000e+03,1.150000000000000000e+03,1.200000000000000000e+03,1.250000000000000000e+03,1.300000000000000000e+03,1.350000000000000000e+03,1.400000000000000000e+03,1.450000000000000000e+03,1.500000000000000000e+03,1.550000000000000000e+03,1.600000000000000000e+03,1.650000000000000000e+03,1.700000000000000000e+03,1.750000000000000000e+03,1.800000000000000000e+03,1.850000000000000000e+03,1.900000000000000000e+03,1.950000000000000000e+03,2.000000000000000000e+03,2.100000000000000000e+03,2.200000000000000000e+03,2.300000000000000000e+03,2.400000000000000000e+03,2.500000000000000000e+03,2.600000000000000000e+03,2.700000000000000000e+03,2.800000000000000000e+03,2.900000000000000000e+03,3.000000000000000000e+03,3.100000000000000000e+03,3.200000000000000000e+03,3.300000000000000000e+03,3.400000000000000000e+03,3.500000000000000000e+03,3.600000000000000000e+03,3.700000000000000000e+03,3.800000000000000000e+03,3.900000000000000000e+03,4.000000000000000000e+03,4.100000000000000000e+03,4.200000000000000000e+03,4.300000000000000000e+03,4.400000000000000000e+03,4.500000000000000000e+03,4.600000000000000000e+03,4.700000000000000000e+03,4.800000000000000000e+03,4.900000000000000000e+03,5.000000000000000000e+03,5.100000000000000000e+03,5.200000000000000000e+03,5.300000000000000000e+03,5.400000000000000000e+03,5.500000000000000000e+03])
zc=zc[::-1]


saltavg=np.array(avgfile.salt.isel(s_rho=lst).isel(eta_rho=ltlst).isel(xi_rho=lglst).groupby('ocean_time.month').mean().values)

a,b=np.meshgrid(lglst, ltlst) 

saltavg=saltavg.mean(axis=0)


for j in range(x_roms.shape[0]):
  for k in range(x_roms.shape[1]):
    if (zr[-1,j,k] > zc.min()):
      zr[-1,j,k] = zc.min()

zr=zr[:,ltmin-ktl:ltmin+ktl, lgmin-ktl:lgmin+ktl]
x_roms=x_roms[ltmin-ktl:ltmin+ktl, lgmin-ktl:lgmin+ktl]
y_roms=y_roms[ltmin-ktl:ltmin+ktl, lgmin-ktl:lgmin+ktl]

isalt=np.zeros([len(zc),x_roms.shape[0], x_roms.shape[1]])


UNDEF=np.nan

for j in range(isalt.shape[1]):
  for k in range(isalt.shape[2]):
    isalt[:,j,k] = np.interp(-zc, -zr[:,j,k], saltavg[:,j,k], right=UNDEF, left=UNDEF)

isalt=np.squeeze(isalt)

saltroms=np.zeros([len(zc)])

for i in range(len(vsitu)):
  saltroms[i]=griddata((x_roms.ravel(),y_roms.ravel()), isalt[i,:].ravel(), (lon,lat))


saltroms=saltroms[::-1]



saltwoa=saltwoa.filled() 
rnan=np.where(isnan(saltroms))
saltwoa[rnan]=np.nan
saltwoa[saltwoa>100]=np.nan


plt.close()
plt.style.use('ggplot')
fig, ax1 = plt.subplots(figsize=(5,7))
ax1.plot(saltroms,-file['depth'][:], linewidth=2,color='blue', label='ROMS')
ax1.plot(saltwoa,-file['depth'][:], linewidth=2, color='red', label='WOA')
ax1.set_xlabel('Salinidade', fontsize=14)
ax1.set_ylabel('Profundidade (m)',fontsize=14)
ax1.grid(alpha=0.5)
legend=ax1.legend(loc=4, fancybox=True, framealpha=1, shadow=True, borderpad=1, fontsize='small')
ax1.set_ylim([-4200,0])
plt.savefig('woa_salt.png', dpi=200, bbox_inches='tight', linestyle='.', transparent=False)
