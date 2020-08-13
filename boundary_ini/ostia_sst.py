cd C:\Users\Fernando\Downloads
activate
python

from netCDF4 import Dataset as dat
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import dates
file=dat('METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2_1590260341312.nc')

lat=np.flipud(file['lat'][:])
lon=file['lon'][:]
time=file['time'][:]
str=np.str(file['time'].units[14:])
num = dates.datestr2num(str) + (time/(24*60*60))
dates.num2date(num)
sst=np.flipud(np.squeeze(file['analysed_sst'][:])) - 273.15
lon,lat=np.meshgrid(lon,lat)



#### trimming

fname_grd = 'azul_grd_era_NEW_menor_azul.nc'  ## Smooth topography grid.   
#fname_grd = 'azul_grd2.nc'  ## Smooth topography grid.   


## Load ROMS grid.
grd = dat(fname_grd)
x_roms = grd.variables['lon_rho'][:]
y_roms = grd.variables['lat_rho'][:]
msk_roms = grd.variables['mask_rho'][:]
msk_romsv = grd.variables['mask_v'][:]

h_roms = grd.variables['h'][:]
sh2 = msk_roms.shape
etamax, ximax = sh2


minlon = lon[0,:] - x_roms.min()
iml = np.where(np.absolute(minlon)==np.absolute(minlon).min())[0][0]
maxlon = lon[0,:] - x_roms.max()
imxl = np.where(np.absolute(maxlon)==np.absolute(maxlon).min())[0][0]

minlat = lat[:,0] - y_roms.min()
imla = np.where(np.absolute(minlat)==np.absolute(minlat).min())[0][0]
maxlat = lat[:,0] - y_roms.max()
imxla = np.where(np.absolute(maxlat)==np.absolute(maxlat).min())[0][0]

lon = lon[imxla:imla,iml:imxl]

lat = lat[imxla:imla,iml:imxl]

sst=sst[imxla:imla,iml:imxl]

plt.pcolor(lon,lat,sst,vmin=20,vmax=27);plt.colorbar();plt.show()