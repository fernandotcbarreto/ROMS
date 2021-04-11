

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

avgfile = xr.open_mfdataset(lista, concat_dim='ocean_time')


hh=np.array(avgfile['ocean_time'][:].dt.strftime('%Y-%m-%d  %H:00:00'))

hhdates=dates.datestr2num(hh)

fname_grd = 'BRSE_2012_GRD.nc'     #CF 1/36

tlst=list(np.arange(0,hh.shape[0],5))

## Load ROMS grid.
grd = Dataset(fname_grd)
x_roms = grd.variables['lon_rho'][:]
y_roms = grd.variables['lat_rho'][:]
msk_roms = grd.variables['mask_rho'][:]
msk_roms=msk_roms.filled()
msk_romsv = grd.variables['mask_v'][:]
h_roms = grd.variables['h'][:]
sh2 = msk_roms.shape
etamax, ximax = sh2


theta_b = 0.4
theta_s = 5.0
tcline = 3.
klevels = 30
Vtransform = 1
Vstretching = 1
Spherical = True

lst=list(np.arange(25,30))

if Vstretching==4:
  scoord = s_coordinate_4(h_roms, theta_b, theta_s, tcline, klevels)   #zeta is not used in the computation 
elif Vstretching==2:
  scoord = s_coordinate_2(h_roms, theta_b, theta_s, tcline, klevels)
elif Vstretching==1:
  scoord = s_coordinate(h_roms, theta_b, theta_s, tcline, klevels)

zr = -scoord.z_r[:]

zr=zr[lst]

zc=np.array([0])



zc=zc[::-1]


tempavg=np.array(avgfile.temp.isel(s_rho=lst).isel(ocean_time=tlst).values)


itemp=np.zeros([tempavg.shape[0],len(zc),x_roms.shape[0], x_roms.shape[1]])

#zeta=avgfile['zeta'][:]

for j in range(itemp.shape[2]):
  for k in range(itemp.shape[3]):
    if (zr[-1,j,k] > zc.min()):
      zr[-1,j,k] = zc.min()


UNDEF=np.nan

for i in range(itemp.shape[0]):
  for j in range(itemp.shape[2]):
    for k in range(itemp.shape[3]):
      itemp[i,:,j,k] = np.interp(-zc, -zr[:,j,k], tempavg[i,:,j,k], right=UNDEF, left=UNDEF)


itemp=np.squeeze(itemp)

import scipy.io as io
#mdic = {"itemp": itemp}
#io.savemat('itemp_diaria.mat', mdic)
itemp=io.loadmat('itemp_diaria.mat')


romstime=hhdates[tlst]

###############MUR
#source: https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41.html


lista = sorted(glob.glob('/mnt/c/Users/Gabriel/Desktop/boundary_ini/jpl_year*'))

murfile=xr.open_mfdataset(lista, concat_dim='time')

murlat=murfile['latitude'][:].values
murlon=murfile['longitude'][:].values

murtemp=np.array(murfile.analysed_sst.values)


hhm=np.array(murfile['time'][:].dt.strftime('%Y-%m-%d  %H:00:00'))

timemur=dates.datestr2num(hhm)


intunewtime=np.zeros([len(timemur),x_roms.shape[0], x_roms.shape[1]])

intvnewtime=np.zeros([len(timemur),x_roms.shape[0], x_roms.shape[1]])

tempnewtime=np.zeros([len(timemur),x_roms.shape[0], x_roms.shape[1]])


for j in range(itemp.shape[1]):
  for k in range(itemp.shape[2]):
#    intunewtime[:,j,k] = np.interp(timemur, romstime, intu[:,j,k], right=UNDEF, left=UNDEF)
#    intvnewtime[:,j,k] = np.interp(timemur, romstime, intv[:,j,k], right=UNDEF, left=UNDEF)
    tempnewtime[:,j,k] = np.interp(timemur, romstime, itemp[:,j,k])


tempinterpmur=np.zeros(tempnewtime.shape)

from scipy.interpolate import interp2d

murtemp=np.ma.masked_invalid(murtemp)
murtemp = np.ma.filled(murtemp,fill_value=100000000)
for i in range(tempinterpmur.shape[0]):
  tempinput = interp2d(murlon, murlat,murtemp[i,:])  
  tempinterpmur[i,:] = tempinput(x_roms[0,:], y_roms[:,0])

#tempnewtime[tempnewtime>100]=np.nan
#tempinterpmur[tempinterpmur>100]=np.nan


#tempdif=tempnewtime-tempinterpmur 

  
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error

temprmsemur=np.zeros(x_roms.shape)
tempmaemur=np.zeros(x_roms.shape)


tempdif=tempnewtime-tempinterpmur 

tempinterpmur=np.ma.masked_invalid(tempinterpmur)
tempnewtime=np.ma.masked_invalid(tempnewtime)
tempinterpmur=tempinterpmur.filled(fill_value=100)
tempnewtime=tempnewtime.filled(fill_value=100)

for j in range(itemp.shape[1]):
  for k in range(itemp.shape[2]):
    temprmsemur[j,k] = np.sqrt(mean_squared_error(tempinterpmur[:,j,k],tempnewtime[:,j,k]))
    tempmaemur[j,k] = mean_absolute_error(tempinterpmur[:,j,k],tempnewtime[:,j,k])

temprmsemur[temprmsemur>5]=5
tempmaemur[tempmaemur>5]=5

#temprmsemur = np.ma.masked_where(temprmsemur>1.5, temprmsemur)
#tempmaemur = np.ma.masked_where(tempmaemur>1.5, tempmaemur)



from mpl_toolkits.basemap import Basemap
import matplotlib.gridspec as gridspec

nlim=y_roms.max()
slim=y_roms.min()
wlim=x_roms.min()
elim=x_roms.max()


temprmsemur[np.where(msk_roms==0)]=np.nan
tempmaemur[np.where(msk_roms==0)]=np.nan
temprmsemur=np.ma.masked_invalid(temprmsemur)
tempmaemur=np.ma.masked_invalid(tempmaemur)


fig = plt.figure(figsize=(10,8))

gs = gridspec.GridSpec(ncols=2, nrows=1)
gs.update(wspace=0.15, hspace=0.1)
ax1 = fig.add_subplot(gs[0,0])
ax3 = fig.add_subplot(gs[0,1])

map = Basemap(projection='merc', llcrnrlat=slim, urcrnrlat=nlim,llcrnrlon=wlim, urcrnrlon=elim,resolution='l', ax=ax1)
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.fillcontinents(color='#ddaa66',lake_color='aqua')
parallels = np.arange(-60,-0,3)
map.drawparallels(parallels,labels=[1,0,0,1], linewidth=0.0, fontsize=7)
meridians = np.arange(-60,-10,4)
map.drawmeridians(meridians,labels=[1,0,0,1], linewidth=0.0, fontsize=7)
#map.readshapefile('/mnt/c/Users/Fernando/Desktop/shape_unid_fed/lim_unidade_federacao_a', 'lim_unidade_federacao_a')
x, y = map(x_roms,y_roms)
a=ax1.pcolor(x,y, np.squeeze(temprmsemur),vmax=3)
#a=ax1.pcolor(x,y, np.squeeze(temprmsemur))
vv,bb=map(-51,-18)
ax1.text(vv,bb,'RMSE\nMean = '+str(np.squeeze(temprmsemur).mean().round(3)), fontsize=9, fontweight='bold', bbox={'facecolor': 'grey', 'alpha': 0.8, 'pad': 3}, horizontalalignment='left',verticalalignment='center')  
vv,bb=map(-55,-15.5)
ax1.text(vv,bb,'A)', fontsize=10, fontweight='bold')

map = Basemap(projection='merc', llcrnrlat=slim, urcrnrlat=nlim,llcrnrlon=wlim, urcrnrlon=elim,resolution='l', ax=ax3)
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.fillcontinents(color='#ddaa66',lake_color='aqua')
parallels = np.arange(-60,-0,3)
map.drawparallels(parallels,labels=[1,0,0,1], linewidth=0.0, fontsize=7)
meridians = np.arange(-60,-10,4)
map.drawmeridians(meridians,labels=[1,0,0,1], linewidth=0.0, fontsize=7)
#map.readshapefile('/mnt/c/Users/Fernando/Desktop/shape_unid_fed/lim_unidade_federacao_a', 'lim_unidade_federacao_a')
x, y = map(x_roms,y_roms)
a=ax3.pcolor(x,y, np.squeeze(tempmaemur),vmax=3)
#a=ax3.pcolor(x,y, np.squeeze(tempmaemur))
cax = fig.add_axes([0.91, 0.3, 0.02, 0.38])
cbar=fig.colorbar(a, shrink=0.6, extend='both', ax=ax3, cax=cax)
cbar.ax.set_ylabel('Temperatura ($^\circ$C)', rotation=270)
cbar.ax.get_yaxis().labelpad = 15
cbar.ax.tick_params(labelsize=8)
vv,bb=map(-51,-18)
ax3.text(vv,bb,'MAE\nMean = '+str(np.squeeze(tempmaemur).mean().round(3)), fontsize=9, fontweight='bold', bbox={'facecolor': 'grey', 'alpha': 0.8, 'pad': 3}, horizontalalignment='left',verticalalignment='center')  
vv,bb=map(-55,-15.5)
ax3.text(vv,bb,'B)', fontsize=10, fontweight='bold')
#plt.savefig('sst_azul', dpi=500) 
plt.savefig('statistics_mur', dpi=300, transparent=False, bbox_inches="tight")
plt.close()

plt.show()



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
cbar.ax.set_ylabel('RMSE Temperature ($^/circ$C)', rotation=270)
cbar.ax.get_yaxis().labelpad = 20
plt.show()



#####################COMPARAÇÃO MENSAL


from scipy.interpolate import interp2d
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
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader
from cartopy.io.shapereader import Reader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import pyPro.utils as pro
import pyPro.maps as mp


lista = sorted(glob.glob('/mnt/share/Modelos/BRSE_2014_2016/RESULTADOS/ocean_BRSE_his_b1_0[0-9][0-9][0-9]*'))
#lista = sorted(glob.glob('/mnt/share/Modelos/BRSE_2014_2016/RESULTADOS/ocean_BRSE_his_b1_00[0-6][0-9]*'))


avgfile = xr.open_mfdataset(lista, concat_dim='ocean_time')
#avgfile = xr.open_mfdataset(lista, concat_dim='ocean_time', decode_cf=False)


hh=np.array(avgfile['ocean_time'][:].dt.strftime('%Y-%m-%d  %H:00:00'))

hhdates=dates.datestr2num(hh)

fname_grd = 'BRSE_2012_GRD.nc'     #CF 1/36

tlst=list(np.arange(0,hh.shape[0],6))

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

zc=np.array([0])



zc=zc[::-1]


tempavg=np.array(avgfile.temp.isel(s_rho=lst).groupby('ocean_time.month').mean().values)


itemp=np.zeros([tempavg.shape[0],len(zc),x_roms.shape[0], x_roms.shape[1]])

#zeta=avgfile['zeta'][:]

for j in range(itemp.shape[2]):
  for k in range(itemp.shape[3]):
    if (zr[-1,j,k] > zc.min()):
      zr[-1,j,k] = zc.min()

UNDEF=np.nan

for i in range(itemp.shape[0]):
  for j in range(itemp.shape[2]):
    for k in range(itemp.shape[3]):
      itemp[i,:,j,k] = np.interp(-zc, -zr[:,j,k], tempavg[i,:,j,k], right=UNDEF, left=UNDEF)


itemp=np.squeeze(itemp)

#import scipy.io as io
#mdic = {"itemp": itemp}
#io.savemat('itemp.mat', mdic)

#MURRR



lista = sorted(glob.glob('/mnt/c/Users/Gabriel/Desktop/boundary_ini/jpl_year*'))

murfile=xr.open_mfdataset(lista, concat_dim='time')

murlat=murfile['latitude'][:].values
murlon=murfile['longitude'][:].values

murtemp=np.array(murfile.analysed_sst.groupby('time.month').mean().values)


tempinterpmur = np.zeros(itemp.shape)
murtemp=np.ma.masked_invalid(murtemp)
murtemp = np.ma.filled(murtemp,fill_value=100000000)
for i in range(tempinterpmur.shape[0]):
  tempinput = interp2d(murlon, murlat,murtemp[i,:], bounds_error=False)  
  tempinterpmur[i,:] = tempinput(x_roms[0,:], y_roms[:,0])

tempinterpmur[tempinterpmur>100]=np.nan


nlim=y_roms.max()
slim=y_roms.min()
wlim=x_roms.min()
elim=x_roms.max()

for ind in range(tempinterpmur.shape[0]):
  if ind == 0:
      tit = 'Janeiro'
  elif ind == 1:
      tit = 'Fevereiro'
  elif ind == 2:
      tit = 'Março'
  elif ind == 3:
      tit = 'Abril'
  elif ind == 4:
      tit = 'Maio'
  elif ind == 5:
      tit = 'Junho'
  elif ind == 6:
      tit = 'Julho'
  elif ind == 7:
      tit = 'Agosto'
  elif ind == 8:
      tit = 'Setembro'
  elif ind == 9:
      tit = 'Outubro'
  elif ind == 10:
      tit = 'Novembro'
  elif ind == 11:
      tit = 'Dezembro'
  fig, ax1 = plt.subplots(2,1, figsize=(10,15), subplot_kw={'projection': ccrs.PlateCarree()})
  fig.subplots_adjust(hspace=0.10, wspace=-0.20)
  crs = ccrs.PlateCarree()
  region = (wlim, elim, slim, nlim) #list = [W, E, S, N]    
  ax1[0].set_extent(region, crs=crs)
  ax1[0].set_title('TSM ($^\circ$C)| ROMS | %s' % (tit),fontsize=14)   
  im = ax1[0].pcolor(x_roms, y_roms, itemp[ind,:],vmin=10,vmax=28, transform=crs)
  cbaxes = fig.add_axes([0.92, 0.3, 0.03, 0.4]) #[x,y,dx,dy]
  cb=plt.colorbar(im, format='%0.2f', cax = cbaxes)
  cb.ax.tick_params(labelsize=13)   
  ax1[0].add_feature(cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor='0.75'))
  ax1[0].coastlines(resolution='10m') # “110m”, “50m”, and “10m”.
  ax1[0].add_feature(cartopy.feature.BORDERS, linestyle='-', linewidths=1, alpha=.5)
  states = cfeature.NaturalEarthFeature(category='cultural', scale='10m', facecolor='none', lw=0.5, name='admin_1_states_provinces_shp')
  ax1[0].add_feature(states, edgecolor='0.5', zorder=2) 
#  ax1[0].text(-40.7, -19.5,'ES',horizontalalignment='left', verticalalignment='center',zorder=3,color='k',fontsize=10, transform=crs)
#  ax1[0].text(-42.7, -22.5,'RJ',horizontalalignment='left', verticalalignment='center',zorder=3,color='k',fontsize=10, transform=crs)
#  ax1[0].text(-43, -19.8,'MG',horizontalalignment='left', verticalalignment='center',zorder=3,color='k',fontsize=10, transform=crs)
  mp.scale_bar(ax1[0], location=[0.425, 0.0425], length=100, linewidth=2) 
  gl = ax1[0].gridlines(draw_labels=True, alpha=0.6, color='0.4', linestyle='dotted')
  gl.xlabels_top = False
  gl.ylabels_left = True
  gl.ylabels_right = False
  gl.xlabel_style = {'size': 12}
  gl.ylabel_style = {'size': 12}
  ax1[1].set_extent(region, crs=crs)
  ax1[1].set_title('TSM ($^\circ$C)| MUR | %s' % (tit),fontsize=16)   
  im = ax1[1].pcolor(x_roms, y_roms, tempinterpmur[ind,:],vmin=10,vmax=28, transform=crs)
  ax1[1].add_feature(cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor='0.75'))
  ax1[1].coastlines(resolution='10m') # “110m”, “50m”, and “10m”.
  ax1[1].add_feature(cartopy.feature.BORDERS, linestyle='-', linewidths=1, alpha=.5)
  states = cfeature.NaturalEarthFeature(category='cultural', scale='10m', facecolor='none', lw=0.5, name='admin_1_states_provinces_shp')
  ax1[1].add_feature(states, edgecolor='0.5', zorder=2) 
#  ax1[1].text(-40.7, -19.5,'ES',horizontalalignment='left', verticalalignment='center',zorder=3,color='k',fontsize=10, transform=crs)
#  ax1[1].text(-42.7, -22.5,'RJ',horizontalalignment='left', verticalalignment='center',zorder=3,color='k',fontsize=10, transform=crs)
#  ax1[1].text(-43, -19.8,'MG',horizontalalignment='left', verticalalignment='center',zorder=3,color='k',fontsize=10, transform=crs)
  mp.scale_bar(ax1[1], location=[0.425, 0.0425], length=100, linewidth=2) 
  gl = ax1[1].gridlines(draw_labels=True, alpha=0.6, color='0.4', linestyle='dotted')
  gl.xlabels_top = False
  gl.ylabels_left = True
  gl.ylabels_right = False
  gl.xlabel_style = {'size': 12}
  gl.ylabel_style = {'size': 12}
  plt.savefig('roms_temp%4.4dm' % ind, dpi=300, transparent=False, bbox_inches="tight")
    
        
for ind in range(tempinterpmur.shape[0]):
  plt.pcolor(x_roms,y_roms,itemp[ind,:],vmin=10,vmax=28);plt.colorbar()
  plt.savefig('temp1'+str(ind)+'.png', dpi=200, bbox_inches='tight', linestyle='.', transparent=False)
  plt.close()
  plt.pcolor(x_roms,y_roms,tempinterpmur[ind,:], vmin=10,vmax=28);plt.colorbar()
  plt.savefig('temp2'+str(ind)+'.png', dpi=200, bbox_inches='tight', linestyle='.', transparent=False)
  plt.close()

ind=0  
plt.pcolor(x_roms,y_roms,itemp[ind,:],vmin=10,vmax=28);plt.colorbar()
plt.savefig('temp1.png', dpi=200, bbox_inches='tight', linestyle='.', transparent=False)
plt.close()
plt.pcolor(x_roms,y_roms,tempinterpmur[ind,:], vmin=10,vmax=28);plt.colorbar()
plt.savefig('temp2.png', dpi=200, bbox_inches='tight', linestyle='.', transparent=False)
plt.close()


plt.pcolor(murtemp[ind,:],vmin=10,vmax=28);plt.colorbar()
plt.savefig('temp2.png', dpi=200, bbox_inches='tight', linestyle='.', transparent=False)
plt.close()