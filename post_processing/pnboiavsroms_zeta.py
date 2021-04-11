import pandas as pd
from matplotlib import dates
import numpy as np
from utils import weim


#gauge location


lat=-20.1813    #stat oil
lon=-38.5535     

#lat=-29.3467   #F norte
#lon=-49.7250       

#lat=-22.1000   #F norte
#lon=-40.0167      

#lat=-28.6033   #sta marta
#lon=-48.8117    

#lat=-25.6234   #galheta
#lon=-48.3176

#lat=-24.0517   #moela
#lon=-46.2683

#lat=-20.3896   #vitoria
#lon=-40.2834

#lat=-20.3896
#lon=-39.7034  


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

lista = sorted(glob.glob('/mnt/share/Modelos/BRSE_2014_2016/RESULTADOS/ocean_BRSE_his_b1_0[0-9][0-9][0-9]*'))

#lista = sorted(glob.glob('/mnt/share/Modelos/BRSE_2012/RESULTADO/ocean_BRSE_his_0[4-9][0-9][0-9]*'))

#lista = sorted(glob.glob('ocean_nne2_his_pyrom_0001*'))


def select_point(ds):
    return ds.drop(['u','v','ubar','vbar'])


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


lglst=np.arange(lgmin-3, lgmin+4)
ltlst=np.arange(ltmin-3, ltmin+4)


zeta=avgfile.zeta.isel(eta_rho=ltlst).isel(xi_rho=lglst).values
x_roms=avgfile.lon_rho.isel(eta_rho=ltlst).isel(xi_rho=lglst).values
y_roms=avgfile.lat_rho.isel(eta_rho=ltlst).isel(xi_rho=lglst).values


izeta=np.squeeze(zeta)


romstime=hhdates[tlst]

vzeta=np.zeros(len(romstime))

for i in range(len(vzeta)):
  vzeta[i]=griddata((x_roms.ravel(),y_roms.ravel()), izeta[i,:].ravel(), (lon,lat))


import scipy.io as io
ff={'zeta':vzeta} 
io.savemat('zeta_stat.mat', ff)

######## tidegaufe
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cbook as cbook
import scipy.io as io

#tidez=pd.read_csv('/mnt/c/Users/Fernando/Downloads/delftdashboard_v02.03.16824/delftdashboard/bin/PORTODEVITORIAIHO.tek', delimiter=',',skiprows=5, header=None)

tidez=pd.read_csv('PORTODEVITORIAIHO.tek', delimiter=',',skiprows=5, header=None)

#tidez=io.loadmat('/mnt/c/Users/Fernando/Downloads/delftdashboard_v02.03.16824/delftdashboard/bin/PORTODEVITORIAIHO.mat')

fig, ax = plt.subplots(figsize=[12,7])
ax.plot(pd.to_datetime(hh[0:tidez.shape[0]], infer_datetime_format=True),tidez[2], 'red', label='Tide gauge', linewidth=1)
ax.plot(pd.to_datetime(hh[0:tidez.shape[0]], infer_datetime_format=True),vzeta[0:tidez.shape[0]],'blue', label='ROMS')
legend=ax.legend(loc=2, fontsize='x-small')
legend.get_frame().set_facecolor('grey')
fmt_half_year = mdates.MonthLocator(interval=1)
ax.xaxis.set_major_locator(fmt_half_year)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%y-%m'))
ax.margins(x=0)

#plt.show()

plt.savefig('zeta_2013.png', bbox_inches='tight');plt.close()
