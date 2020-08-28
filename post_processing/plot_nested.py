from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as plt
from netCDF4 import Dataset 
import numpy as np
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.patches as patches
import matplotlib.patheffects as PathEffects
import cmocean


fname_grd = 'grid_rotated_SUL_2.nc'
grd = Dataset(fname_grd)
x_roms = grd.variables['lon_rho'][:]
y_roms = grd.variables['lat_rho'][:]
h_roms = grd.variables['h'][:]

fname_grd = 'grid_rotated_SUL_2_NEST_smaler.nc'
grd2 = Dataset(fname_grd)
x_roms2 = grd2.variables['lon_rho'][:]
y_roms2 = grd2.variables['lat_rho'][:]


##great bath
#batpath=Dataset('/mnt/c/Users/Fernando/Desktop/make_grid/gebco_2020.nc')
#latbat=np.ma.filled(batpath['lat'][:])
#lonbat=np.ma.filled(batpath['lon'][:])
#batgt=np.ma.filled(batpath['elevation'][:])
#lonbat,latbat=np.meshgrid(lonbat,latbat)

batpath=Dataset('/mnt/c/Users/Fernando/Desktop/make_grid/ETOPO2v2g_f4.nc')
latbat=np.ma.filled(batpath['y'][:])
lonbat=np.ma.filled(batpath['x'][:])
batgt=np.ma.filled(batpath['z'][:])
lonbat,latbat=np.meshgrid(lonbat,latbat)
batgt=batgt*-1
batgt[batgt<0]=np.nan
batgt=np.ma.masked_invalid(batgt)

row,col=np.where((latbat<-10) & (latbat>-40) & (lonbat>-60) & (lonbat<-20))
row=np.unique(row)
col=np.unique(col)
col,row=np.meshgrid(col,row)
latbat=latbat[row,col]
lonbat=lonbat[row,col]
batgt=batgt[row,col]
#plt.pcolor(lonbat,latbat, batgt, cmap=cmocean.cm.cmap_d['deep']);plt.show()


lon=x_roms.copy()
lat=y_roms.copy()

lon2=x_roms2.copy()
lat2=y_roms2.copy()

nlim=-20.7+4
slim=lat.min()-5
wlim=lon.min()-5
elim=lon.max()+7




fig, ax = plt.subplots(figsize=(15,15))
map = Basemap(projection='merc', llcrnrlat=slim, urcrnrlat=nlim,llcrnrlon=wlim, urcrnrlon=elim,resolution='l')
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.fillcontinents(color='silver',lake_color='white')
parallels = np.arange(-40,-15,3)
map.drawparallels(parallels,labels=[1,0,0,1], linewidth=0.4)
meridians = np.arange(-55,-30,3)
map.drawmeridians(meridians,labels=[1,0,0,1], linewidth=0.4)
map.readshapefile('/mnt/c/Users/Fernando/Desktop/shape_unid_fed/lim_unidade_federacao_a', 'lim_unidade_federacao_a')
map.drawmapscale(-37., -37, -37, -37, 500, barstyle='fancy', fontsize = 8, yoffset=20000)
x, y = map(lon,lat)
#ax.plot(x,y, 'blue')
x2, y2 = map(lon2,lat2)
#ax.plot(x2,y2, 'red')
x1,y1 = map(lon[0,:].min()+0.5,lat[0,:].max()+0.5)
x2,y2 = map(lon[-1,:].min()+0.5,lat[-1,:].max()+0.5)
x3,y3 = map(lon[-1,:].max(),lat[-1,:].min()+0.5)
x4,y4 = map(lon[0,:].max(),lat[0,:].min()+0.5)
poly = Polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4)],facecolor='none',edgecolor='navy',linewidth=2, alpha=1)
ax.add_patch(poly)
x1,y1 = map(-55,-34)
x2,y2 = map(-48.5,-18)
x3,y3 = map(-35,-22.5)
x4,y4 = map(-41,-38)
poly = Polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4)],facecolor='none',edgecolor='darkred',linewidth=2, alpha=1)
ax.add_patch(poly)
#cnt=ax.contour(x,y,h_roms,levels=[200, 2000],colors=('k','red'))
x,y=map(lonbat,latbat)
cnt=ax.contour(x,y,batgt,levels=[100,2000,3500],colors=('k'), linewidths=1)
ax.pcolor(x,y,batgt,cmap=cmocean.cm.cmap_d['deep'])
clabels=ax.clabel(cnt, inline=True, colors='k', fmt = '%2.1d',fontsize=7)
[txt.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')]) for txt in clabels]
x,y=map(-47.15,-27.24)
ax.plot(x,y,color='indigo', linestyle='none', marker="X", markersize=9, label='Boia Itajaí')
ax.legend(loc=1, bbox_to_anchor=(0.25,0.1), prop={'weight':'bold'})
for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(1)

#plt.show()

#left, bottom, width, height = [0.6, 0.6, 0.2, 0.2]
#ax2 = fig.add_axes([left, bottom, width, height])

#map = Basemap(projection='merc', llcrnrlat=slim, urcrnrlat=nlim,llcrnrlon=wlim, urcrnrlon=elim,resolution='l')
#map.drawcoastlines(linewidth=0.25)
#map.drawcountries(linewidth=0.25)
#map.fillcontinents(color='#ddaa66',lake_color='aqua')
#parallels = np.arange(-40,-15,2)
#map.drawparallels(parallels,labels=[1,0,0,1], linewidth=0.0)
#meridians = np.arange(-55,-30,2)
#map.drawmeridians(meridians,labels=[1,0,0,1], linewidth=0.0)
#map.readshapefile('/mnt/c/Users/Fernando/Desktop/shape_unid_fed/lim_unidade_federacao_a', 'lim_unidade_federacao_a')
#x, y = map(lon,lat)
#ax2.plot(x,y, 'blue')


left, bottom, width, height = [0.25, 0.8, 0.15, 0.15]
ax3 = fig.add_axes([left, bottom, width, height])

map = Basemap(projection='ortho',lat_0=-20,lon_0=-50,resolution='l', ax=ax3)
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.fillcontinents(color='grey',lake_color='aqua')
map.drawmapboundary(fill_color='powderblue')
#map.drawmeridians(np.arange(0,360,30))
#map.drawparallels(np.arange(-90,90,30))
#x, y = map(lon,lat)
#ax2.plot(x,y, 'blue')


x1,y1 = map(-55,-35)
x2,y2 = map(-55,-19)
x3,y3 = map(-35,-19)
x4,y4 = map(-35,-35)
poly = Polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4)],facecolor='blue',edgecolor='black',linewidth=1, alpha=0.8)
ax3.add_patch(poly)

plt.tight_layout()

plt.show()

plt.gca().set_axis_off()
plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, 
            hspace = 0, wspace = 0)
plt.margins(0,0)
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.savefig("itajai.png", bbox_inches = 'tight',pad_inches = 0)

plt.show()

