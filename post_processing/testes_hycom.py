
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from vgrid import s_coordinate, s_coordinate_2, s_coordinate_4
import numpy.ma as ma
from scipy.interpolate import griddata, interp1d
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as dat
import time
import pandas as pd
from matplotlib.pylab import *
from timecalc import timeit
import matplotlib.dates as dates
from mpl_toolkits.basemap import Basemap






###############análise myocean

nc_hycom='/home/fernando/roms/src/Projects/operational/mercator_data/MYOCEAN_AZUL_FORECAST_20200816.nc'

##nc_hycom='/home/fernando/roms/src/Projects/hindcast_2/mercator/MYOCEAN_AZUL_FORECAST_20171124.nc'

file=Dataset(nc_hycom)

x_fm=file['longitude'][:]

y_fm=file['latitude'][:]
y_fm=y_fm[::-1] 

x_fm,y_fm=np.meshgrid(x_fm,y_fm)

z_fm = np.flipud(-file['depth'][:])   

temp_fm=np.squeeze(file['thetao'][:])
temp_fm=temp_fm.filled()
temp_fm[temp_fm<-100]=np.nan
temp_fm=temp_fm[::-1,::-1,:]

salt_fm=np.squeeze(file['so'][:])
salt_fm=salt_fm.filled()
salt_fm[salt_fm<-100]=np.nan
salt_fm=salt_fm[::-1,::-1,:]

u_fm=np.squeeze(file['uo'][:])
u_fm=u_fm.filled()  
u_fm[u_fm<-100]=np.nan
u_fm=u_fm[::-1,::-1,:]

v_fm=np.squeeze(file['vo'][:])
v_fm=v_fm.filled()  
v_fm[v_fm<-100]=np.nan  
v_fm=v_fm[::-1,::-1,:]

ssh_fm = np.squeeze((file['zos'][:]))
ssh_fm=ssh_fm.filled()  
ssh_fm[ssh_fm<-100]=np.nan  
ssh_fm=ssh_fm[::-1,:]

fname_grd = 'grid_rotated_SUL_2_smaler.nc'
#fname_grd = 'azul_grd2.nc'  ## Smooth topography grid.   


## Load ROMS grid.
grd = Dataset(fname_grd)
x_roms = grd.variables['lon_rho'][:]
y_roms = grd.variables['lat_rho'][:]
msk_roms = grd.variables['mask_rho'][:]
msk_romsv = grd.variables['mask_v'][:]

h_roms = grd.variables['h'][:]
sh2 = msk_roms.shape
etamax, ximax = sh2


minlon = x_fm[0,:] - x_roms.min()
iml = np.where(np.absolute(minlon)==np.absolute(minlon).min())[0][0]
maxlon = x_fm[0,:] - x_roms.max()
imxl = np.where(np.absolute(maxlon)==np.absolute(maxlon).min())[0][0]

minlat = y_fm[:,0] - y_roms.min()
imla = np.where(np.absolute(minlat)==np.absolute(minlat).min())[0][0]
maxlat = y_fm[:,0] - y_roms.max()
imxla = np.where(np.absolute(maxlat)==np.absolute(maxlat).min())[0][0]

x_fm = x_fm[imxla:imla,iml:imxl]

y_fm = y_fm[imxla:imla,iml:imxl]

temp_fm = temp_fm[:,imxla:imla,iml:imxl]
salt_fm = salt_fm[:,imxla:imla,iml:imxl]
u_fm = u_fm[:,imxla:imla,iml:imxl]
v_fm = v_fm[:,imxla:imla,iml:imxl]
ssh_fm = ssh_fm[imxla:imla,iml:imxl]


###############################



lay=28

sp=2


uc=ma.masked_invalid(np.squeeze(u_fm[lay,:,:]))
vc=ma.masked_invalid(np.squeeze(v_fm[lay,:,:]))
val=np.sqrt((uc**2)+(vc**2))

a=plt.pcolor(x_fm, y_fm,val, vmin=val.min(), vmax=0.8)
#a=plt.pcolor(x_fm, y_fm,val, vmin=val.min(), vmax=val.max())
plt.quiver(x_fm[0:-1:sp, 0:-1:sp],y_fm[0:-1:sp, 0:-1:sp],uc[0:-1:sp, 0:-1:sp], vc[0:-1:sp, 0:-1:sp])
plt.colorbar(a)
plt.show()




plt.pcolor(x_fm,y_fm, vc, vmax=0.5, vmin=-1.3);plt.colorbar();plt.show()

plt.pcolor(x_fm,y_fm, uc, vmax=0.5, vmin=-0.7);plt.colorbar();plt.show()



lay=-1

a=ma.masked_invalid(temp_fm[lay,::])

plt.pcolor(x_fm, y_fm,a, vmin=15,vmax=24)
plt.colorbar()
plt.show()





# NORTH SOUTH

a,b=np.meshgrid(x_fm[0,:],z_fm)

fig, ax1 = plt.subplots()
AB=ax1.pcolor(a,b,ma.masked_invalid(v_fm[:,-2,:]), vmin=-1, vmax=1, cmap = 'seismic')
ax1.axis('tight')
cbar=fig.colorbar(AB, shrink=0.9, extend='both')
cbar.ax.set_ylabel('v($m^{2}$)', rotation=270)
cbar.ax.get_yaxis().labelpad = 20
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=15)
text.set_font_properties(font)
#plt.savefig('myocean2.png')
plt.ylim(-5000,0)
plt.show()



######################   EAST WEST

a,b=np.meshgrid(y_fm[:,-1],z_fm)

fig, ax1 = plt.subplots()
AB=ax1.pcolor(a,b,ma.masked_invalid(u_fm[:,:,-3]), vmin=-0.15, vmax=0.15)
ax1.axis('tight')
cbar=fig.colorbar(AB, shrink=0.9, extend='both')
cbar.ax.set_ylabel('v($m^{2}$)', rotation=270)
cbar.ax.get_yaxis().labelpad = 20
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=15)
text.set_font_properties(font)
#plt.savefig('myocean2.png')

plt.show()









file=dat('prooceano_myocean_cortado_Forecast_jj_bry.nc')
temp=file['v_south'][::]

theta_b = 0.4
theta_s = 5.0
tcline = 100.
klevels = 40
Vtransform = 2
Vstretching = 4
fname_grd='azul_grd2.nc'
grd = dat(fname_grd)
x_roms = grd.variables['lon_rho'][:]
y_roms = grd.variables['lat_rho'][:]
msk_roms = grd.variables['mask_rho'][:]
msk_romsv = grd.variables['mask_v'][:]
h_roms = grd.variables['h'][:]
sh2 = msk_roms.shape
etamax, ximax = sh2

scoord = s_coordinate_4(h_roms, theta_b, theta_s, tcline, klevels)
zr = scoord.z_r[:]







#EAST WEST  subplot myocean roms

a,b=np.meshgrid(y_fm[:,-1],z_fm)

fig, ax1 = plt.subplots(1,2, figsize=(10,5))
AB=ax1[0].pcolor(a,b,ma.masked_invalid(v_fm[:,:,-3]), vmin=-0.15, vmax=0.15, cmap = 'RdBu')
#AB=ax1[0].pcolor(a,b,ma.masked_invalid(salt_fm[:,:,-3]))
ax1[0].axis('tight')
cbar=fig.colorbar(AB, shrink=0.9, extend='both', ax=ax1[0])
cbar.ax.set_ylabel('V-velocity(m/s)', rotation=270)
cbar.ax.get_yaxis().labelpad = 20
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=15)
text.set_font_properties(font)
ax1[0].title.set_text('MyOcean')
#plt.savefig('myocean2.png')
ax1[0].set_ylim([-5000,0])

a,b=np.meshgrid(y_roms[:,0],zr[:,0,0])

AB=ax1[1].pcolor(a, zr[:,:,-1], temp[-1,::], vmin=-0.15, vmax=0.15, cmap = 'RdBu')
#AB=ax1[1].pcolor(a, zr[:,:,-1], temp[5,::])
ax1[1].axis('tight')
cbar=fig.colorbar(AB, shrink=0.9, extend='both', ax=ax1[1])
cbar.ax.set_ylabel('V-velocity(m/s)', rotation=270)
cbar.ax.get_yaxis().labelpad = 20
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=15)
ax1[1].title.set_text('Roms')
text.set_font_properties(font)
ax1[1].set_ylim([-5000,0])


plt.tight_layout()

plt.savefig('v_east.png')
plt.show()






#SOUTH NORTH subplot myocean roms



a,b=np.meshgrid(x_fm[0,:],z_fm)

fig, ax1 = plt.subplots(1,2, figsize=(10,5))
AB=ax1[0].pcolor(a,b,ma.masked_invalid(v_fm[:,-3,:]), vmin=-0.3, vmax=0.3, cmap = 'RdBu')
#AB=ax1[0].pcolor(a,b,ma.masked_invalid(temp_fm[:,-3,:]))
ax1[0].axis('tight')
cbar=fig.colorbar(AB, shrink=0.9, extend='both', ax=ax1[0])
cbar.ax.set_ylabel('V-velocity(m/s)', rotation=270)
cbar.ax.get_yaxis().labelpad = 20
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=15)
text.set_font_properties(font)
ax1[0].title.set_text('MyOcean')
#plt.savefig('myocean2.png')
ax1[0].set_ylim([-5000,0])

a,b=np.meshgrid(x_roms[0,:],zr[:,0,0])

AB=ax1[1].pcolor(a, zr[:,0,:], temp[-1,::], vmin=-0.3, vmax=0.3, cmap = 'RdBu')
#AB=ax1[1].pcolor(a, zr[:,0,:], temp[5,::])
ax1[1].axis('tight')
cbar=fig.colorbar(AB, shrink=0.9, extend='both', ax=ax1[1])
cbar.ax.set_ylabel('V-velocity(m/s)', rotation=270)
cbar.ax.get_yaxis().labelpad = 20
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=15)
ax1[1].title.set_text('Roms')
text.set_font_properties(font)
plt.ylim(-5000,0)
ax1[1].set_ylim([-5000,0])


plt.tight_layout()

#plt.savefig('v_south.png')
plt.show()


##############################################subplot myocean roms BASEMAP


nc_hycom='C:\Users\Fernando\Desktop\dados_myocean\MYOCEAN_AZUL_FORECAST_20200423.nc'

file=Dataset(nc_hycom)

x_fm=file['longitude'][:]

y_fm=file['latitude'][:]
y_fm=y_fm[::-1] 

x_fm,y_fm=np.meshgrid(x_fm,y_fm)

z_fm = np.flipud(-file['depth'][:])   

temp_fm=np.squeeze(file['thetao'][:])
temp_fm=temp_fm.filled()
temp_fm[temp_fm<-100]=np.nan
temp_fm=temp_fm[::-1,::-1,:]

salt_fm=np.squeeze(file['so'][:])
salt_fm=salt_fm.filled()
salt_fm[salt_fm<-100]=np.nan
salt_fm=salt_fm[::-1,::-1,:]

u_fm=np.squeeze(file['uo'][:])
u_fm=u_fm.filled()  
u_fm[u_fm<-100]=np.nan
u_fm=u_fm[::-1,::-1,:]

v_fm=np.squeeze(file['vo'][:])
v_fm=v_fm.filled()  
v_fm[v_fm<-100]=np.nan  
v_fm=v_fm[::-1,::-1,:]

ssh_fm = np.squeeze((file['zos'][:]))
ssh_fm=ssh_fm.filled()  
ssh_fm[ssh_fm<-100]=np.nan  
ssh_fm=ssh_fm[::-1,:]

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


minlon = x_fm[0,:] - x_roms.min()
iml = np.where(np.absolute(minlon)==np.absolute(minlon).min())[0][0]
maxlon = x_fm[0,:] - x_roms.max()
imxl = np.where(np.absolute(maxlon)==np.absolute(maxlon).min())[0][0]

minlat = y_fm[:,0] - y_roms.min()
imla = np.where(np.absolute(minlat)==np.absolute(minlat).min())[0][0]
maxlat = y_fm[:,0] - y_roms.max()
imxla = np.where(np.absolute(maxlat)==np.absolute(maxlat).min())[0][0]

x_fm = x_fm[imxla-2:imla+2,iml-2:imxl+2]

y_fm = y_fm[imxla-2:imla+2,iml-2:imxl+2]

temp_fm = temp_fm[:,imxla-2:imla+2,iml-2:imxl+2]
salt_fm = salt_fm[:,imxla-2:imla+2,iml-2:imxl+2]
u_fm = u_fm[:,imxla-2:imla+2,iml-2:imxl+2]
v_fm = v_fm[:,imxla-2:imla+2,iml-2:imxl+2]
ssh_fm = ssh_fm[imxla-2:imla+2,iml-2:imxl+2]




sp=4
tim=-1
lay=-1

rj=(-42.7,-22.5)

nlim=lat.max()
slim=lat.min()
wlim=lon.min()
elim=lon.max()



########horizontal plots

fig, ax1 = plt.subplots(1,2, figsize=(12,8))

val=np.sqrt((u[tim,lay,:,:]**2)+(v[tim,lay,:,:]**2))

map1 = Basemap(projection='merc', llcrnrlat=slim, urcrnrlat=nlim,llcrnrlon=wlim, urcrnrlon=elim,resolution='l', ax=ax1[0])
# draw coastlines, country boundaries, fill continents.
map1.drawcoastlines(linewidth=0.25)
map1.drawcountries(linewidth=0.25)
map1.fillcontinents(color='#ddaa66',lake_color='aqua')
parallels = np.arange(-30,-10,2)
map1.drawparallels(parallels,labels=[1,0,0,1], linewidth=0.0)
meridians = np.arange(-55,-30,2)
map1.drawmeridians(meridians,labels=[1,0,0,1], linewidth=0.0)
map1.readshapefile('/mnt/c/Users/Fernando/Desktop/shape_unid_fed/lim_unidade_federacao_a', 'lim_unidade_federacao_a')
x, y = map1(lon,lat)
AB=ax1[0].pcolor(x,y, np.squeeze(val),vmin=val.min(),vmax=val.max())
ax1[0].quiver(x[0:-1:sp, 0:-1:sp],y[0:-1:sp, 0:-1:sp],u[tim,lay,0:-1:sp, 0:-1:sp], v[tim,lay,0:-1:sp, 0:-1:sp])
ax1[0].set_title('Roms')
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=15)
#ax1[0].axis('tight')  #destroi basemap


uc=ma.masked_invalid(np.squeeze(u_fm[lay,:,:]))
vc=ma.masked_invalid(np.squeeze(v_fm[lay,:,:]))
val=np.sqrt((uc**2)+(vc**2))

ax1[1].title.set_text('MyOcean')

map = Basemap(projection='merc', llcrnrlat=slim, urcrnrlat=nlim,llcrnrlon=wlim, urcrnrlon=elim,resolution='l', ax=ax1[1])
# draw coastlines, country boundaries, fill continents.
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.fillcontinents(color='#ddaa66',lake_color='aqua')
parallels = np.arange(-30,-10,2)
map.drawparallels(parallels,labels=[1,0,0,1], linewidth=0.0)
meridians = np.arange(-55,-30,2)
map.drawmeridians(meridians,labels=[1,0,0,1], linewidth=0.0)
map.readshapefile('/mnt/c/Users/Fernando/Desktop/shape_unid_fed/lim_unidade_federacao_a', 'lim_unidade_federacao_a')
x, y = map(x_fm,y_fm)
AB=ax1[1].pcolor(x,y, val,vmin=val.min(), vmax=val.max())
ax1[1].quiver(x[0:-1:sp, 0:-1:sp],y[0:-1:sp, 0:-1:sp],uc[0:-1:sp, 0:-1:sp], vc[0:-1:sp, 0:-1:sp])
#ax1[1].axis('tight')  #destroi basemap
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=15)


p0 = ax1[0].get_position().get_points().flatten()
p1 = ax1[1].get_position().get_points().flatten()

ax_cbar = fig.add_axes([0.8, 0.3, 0.03, 0.4])

cbar = plt.colorbar(AB, cax=ax_cbar)

cbar.ax.set_ylabel('Velocity(m/s)', rotation=270)
cbar.ax.get_yaxis().labelpad = 15
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=10)
text.set_font_properties(font)

#fig.tight_layout()


plt.show()


######## vertical plots

sp1=6
sp=3
tim=-1
lay=35

velmin=0
velmax=0.9

rj=(-42.7,-22.5)

nlim=lat.max()
slim=lat.min()
wlim=lon.min()
elim=lon.max()


fig, ax1 = plt.subplots(2,1, figsize=(8,15))

val=np.sqrt((u[tim,lay,:,:]**2)+(v[tim,lay,:,:]**2))

map1 = Basemap(projection='merc', llcrnrlat=slim, urcrnrlat=nlim,llcrnrlon=wlim, urcrnrlon=elim,resolution='l', ax=ax1[0])
# draw coastlines, country boundaries, fill continents.
map1.drawcoastlines(linewidth=0.25)
map1.drawcountries(linewidth=0.25)
map1.fillcontinents(color='#ddaa66',lake_color='aqua')
parallels = np.arange(-30,-10,2)
map1.drawparallels(parallels,labels=[1,0,0,1], linewidth=0.0)
meridians = np.arange(-55,-25,4)
map1.drawmeridians(meridians,labels=[1,0,0,1], linewidth=0.0)
map1.readshapefile('C:\Users\Fernando\Desktop\shape_unid_fed\lim_unidade_federacao_a', 'lim_unidade_federacao_a')
x, y = map1(lon,lat)
AB=ax1[0].pcolor(x,y, np.squeeze(val),vmin=velmin,vmax=velmax)
ax1[0].quiver(x[0:-1:sp1, 0:-1:sp1],y[0:-1:sp1, 0:-1:sp1],u[tim,lay,0:-1:sp1, 0:-1:sp1], v[tim,lay,0:-1:sp1, 0:-1:sp1])
ax1[0].set_title('Nested Roms 0.02$^\circ$')
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=15)


uc=ma.masked_invalid(np.squeeze(u_fm[lay,:,:]))
vc=ma.masked_invalid(np.squeeze(v_fm[lay,:,:]))
val=np.sqrt((uc**2)+(vc**2))

ax1[1].title.set_text('MyOcean')

map = Basemap(projection='merc', llcrnrlat=slim, urcrnrlat=nlim,llcrnrlon=wlim, urcrnrlon=elim,resolution='l', ax=ax1[1])
# draw coastlines, country boundaries, fill continents.
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.fillcontinents(color='#ddaa66',lake_color='aqua')
parallels = np.arange(-30,-10,2)
map.drawparallels(parallels,labels=[1,0,0,1], linewidth=0.0)
meridians = np.arange(-55,-25,4)
map.drawmeridians(meridians,labels=[1,0,0,1], linewidth=0.0)
map1.readshapefile('C:\Users\Fernando\Desktop\shape_unid_fed\lim_unidade_federacao_a', 'lim_unidade_federacao_a')
x, y = map(x_fm,y_fm)
AB=ax1[1].pcolor(x,y, val,vmin=velmin,vmax=velmax)
ax1[1].quiver(x[0:-1:sp, 0:-1:sp],y[0:-1:sp, 0:-1:sp],uc[0:-1:sp, 0:-1:sp], vc[0:-1:sp, 0:-1:sp])
#ax1[1].axis('tight')  #destroi basemap
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=15)


p0 = ax1[0].get_position().get_points().flatten()
p1 = ax1[1].get_position().get_points().flatten()

ax_cbar = fig.add_axes([0.8, 0.3, 0.03, 0.4])

cbar = plt.colorbar(AB, cax=ax_cbar)

cbar.ax.set_ylabel('Velocity(m/s)', rotation=270)
cbar.ax.get_yaxis().labelpad = 15
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=10)
text.set_font_properties(font)

plt.suptitle('Depth = ' + str (zc[lay]) + ' m')

plt.show()


#TEMPERATURA

tim=-2
lay=-1

rj=(-42.7,-22.5)

nlim=lat.max()
slim=lat.min()
wlim=lon.min()
elim=lon.max()


fig, ax1 = plt.subplots(2,1, figsize=(8,15))


map1 = Basemap(projection='merc', llcrnrlat=slim, urcrnrlat=nlim,llcrnrlon=wlim, urcrnrlon=elim,resolution='l', ax=ax1[0])
# draw coastlines, country boundaries, fill continents.
map1.drawcoastlines(linewidth=0.25)
map1.drawcountries(linewidth=0.25)
map1.fillcontinents(color='#ddaa66',lake_color='aqua')
parallels = np.arange(-30,-10,2)
map1.drawparallels(parallels,labels=[1,0,0,1], linewidth=0.0)
meridians = np.arange(-55,-25,4)
map1.drawmeridians(meridians,labels=[1,0,0,1], linewidth=0.0)
map1.readshapefile('C:\Users\Fernando\Desktop\shape_unid_fed\lim_unidade_federacao_a', 'lim_unidade_federacao_a')
x, y = map1(lon,lat)
AB=ax1[0].pcolor(x,y, np.squeeze(temp[tim,lay,:,:]), vmin=21,vmax=29)
ax1[0].set_title('Roms')
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=15)


uc=ma.masked_invalid(np.squeeze(u_fm[lay,:,:]))
vc=ma.masked_invalid(np.squeeze(v_fm[lay,:,:]))
val=np.sqrt((uc**2)+(vc**2))

ax1[1].title.set_text('MyOcean')

map = Basemap(projection='merc', llcrnrlat=slim, urcrnrlat=nlim,llcrnrlon=wlim, urcrnrlon=elim,resolution='l', ax=ax1[1])
# draw coastlines, country boundaries, fill continents.
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.fillcontinents(color='#ddaa66',lake_color='aqua')
parallels = np.arange(-30,-10,2)
map.drawparallels(parallels,labels=[1,0,0,1], linewidth=0.0)
meridians = np.arange(-55,-25,4)
map.drawmeridians(meridians,labels=[1,0,0,1], linewidth=0.0)
map1.readshapefile('C:\Users\Fernando\Desktop\shape_unid_fed\lim_unidade_federacao_a', 'lim_unidade_federacao_a')
x, y = map(x_fm,y_fm)
AB=ax1[1].pcolor(x,y, temp_fm[lay,:,:],vmin=21,vmax=29)
#ax1[1].axis('tight')  #destroi basemap
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=15)

p0 = ax1[0].get_position().get_points().flatten()
p1 = ax1[1].get_position().get_points().flatten()

ax_cbar = fig.add_axes([0.8, 0.3, 0.03, 0.4])

cbar = plt.colorbar(AB, cax=ax_cbar)

cbar.ax.set_ylabel('Temperatura ($^\circ$C)', rotation=270)
cbar.ax.get_yaxis().labelpad = 15
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=12)
text.set_font_properties(font)

plt.suptitle('Depth = ' + str (zc[lay]) + ' m')

plt.show()




AB=ax1[1].pcolor(a, zr[:,0,:], temp[-1,::], vmin=-0.3, vmax=0.3, cmap = 'RdBu')
#AB=ax1[1].pcolor(a, zr[:,0,:], temp[5,::])
ax1[1].axis('tight')
cbar=fig.colorbar(AB, shrink=0.9, extend='both', ax=ax1[1])
cbar.ax.set_ylabel('V-velocity(m/s)', rotation=270)
cbar.ax.get_yaxis().labelpad = 20
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=15)
ax1[1].title.set_text('Roms')
text.set_font_properties(font)
plt.ylim(-5000,0)
ax1[1].set_ylim([-5000,0])















PLOT BRY

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset as dat
from vgrid import s_coordinate, s_coordinate_2, s_coordinate_4
from matplotlib.pylab import *
from matplotlib import dates



file=dat('prooceano_rotate_azul_bry.nc')


#plt.pcolor(temp[0,:,:]);plt.colorbar();plt.show()

#plt.pcolor(temp1[7,:,:]);plt.colorbar();plt.show()

ddd=dates.num2date(dates.datestr2num('2013-01-01')+file['temp_time'][:]/(24*60*60))

theta_b = 0.4
theta_s = 5.0
tcline = 100.
klevels = 40
Vtransform = 2
Vstretching = 4

fname_grd='grid_rotated_azul.nc'

grd = dat(fname_grd)
x_roms = grd.variables['lon_rho'][:]
y_roms = grd.variables['lat_rho'][:]
msk_roms = grd.variables['mask_rho'][:]
msk_romsv = grd.variables['mask_v'][:]
h_roms = grd.variables['h'][:]
sh2 = msk_roms.shape
etamax, ximax = sh2

scoord = s_coordinate_4(h_roms, theta_b, theta_s, tcline, klevels)
zr = scoord.z_r[:]

a,b=np.meshgrid(x_roms[-1,:],zr[:,-1,0])

temp=file['u_north'][::]

valbry=temp[-1,::]

fig, ax1 = plt.subplots()
AB=ax1.pcolor(a, zr[:,-1,:],valbry,  vmin=valbry.min(), vmax=valbry.max(), cmap = 'seismic')
#AB=ax1.pcolor(a, zr[:,0,:], temp[-1,::],  cmap = 'seismic')
#AB=ax1.pcolor(a, zr[:,-1,:], temp[0,::])
ax1.axis('tight')
cbar=fig.colorbar(AB, shrink=0.9, extend='both')
cbar.ax.set_ylabel('v($m^{2}$)', rotation=270)
cbar.ax.get_yaxis().labelpad = 20
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=15)
text.set_font_properties(font)
#plt.savefig('bry_south.png')
plt.ylim(-5000,0)
plt.show()





plt.show()



valbry=temp[0,::]

fig, ax1 = plt.subplots()
AB=ax1.pcolor(a, zr[:,0,:],valbry,  vmin=valbry.min(), vmax=valbry.max(), cmap = 'seismic')
#AB=ax1.pcolor(a, zr[:,0,:], temp[-1,::],  cmap = 'seismic')
#AB=ax1.pcolor(a, zr[:,-1,:], temp[0,::])
ax1.axis('tight')
cbar=fig.colorbar(AB, shrink=0.9, extend='both')
cbar.ax.set_ylabel('v($m^{2}$)', rotation=270)
cbar.ax.get_yaxis().labelpad = 20
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=15)
text.set_font_properties(font)
#plt.savefig('bry_south.png')
plt.ylim(-5000,0)

(zind,lonind)=np.where(valbry==valbry.max())
ax1.plot(a[zind,lonind], zr[zind,0,lonind], '*', color='green')
(zind,lonind)=np.where(valbry==valbry.min())
ax1.plot(a[zind,lonind], zr[zind,0,lonind], '*', color='black')

plt.show()


########################################2D PLOTS

lay=-1

sp=3


uc=ma.masked_invalid(np.squeeze(u_fm[lay,:,:]))
vc=ma.masked_invalid(np.squeeze(v_fm[lay,:,:]))
val=np.sqrt((uc**2)+(vc**2))

a=plt.pcolor(x_fm, y_fm,val, vmin=val.min(), vmax=1)
plt.quiver(x_fm[0:-1:sp, 0:-1:sp],y_fm[0:-1:sp, 0:-1:sp],uc[0:-1:sp, 0:-1:sp], vc[0:-1:sp, 0:-1:sp])
plt.colorbar(a)
plt.show()



lay=-1

a=ma.masked_invalid(temp_fm[lay,::])

plt.pcolor(x_fm, y_fm,a, vmin=a.min(),vmax=29)
plt.colorbar()
plt.show()