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
from matplotlib.collections import PatchCollection
import matplotlib.patches as patches
import matplotlib.patheffects as PathEffects

#lista = sorted(glob.glob('R:/Modelos/BRSE_2014_2016/RESULTADOS/ocean_BRSE_his_b22_00[2-9][0-9]*')) + sorted(glob.glob('R:/Modelos/BRSE_2014_2016/RESULTADOS/ocean_BRSE_his_b3_00[0-4][0-9]*'))


#lista = sorted(glob.glob('R:/Modelos/BRSE_2014_2016/RESULTADOS/ocean_BRSE_his_b1_0[0-9][0-9][0-9]*'))

lista = sorted(glob.glob('/mnt/share/Modelos/BRSE_2014_2016/RESULTADOS/ocean_BRSE_his_b1_0[0-9][0-9][0-9]*'))

avgfile = xr.open_mfdataset(lista, concat_dim='ocean_time')

hh=np.array(avgfile['ocean_time'][:].dt.strftime('%Y-%m-%d  %H:00:00'))

hhdates=dates.datestr2num(hh)

fname_grd = 'BRSE_2012_GRD.nc'     #CF 1/36

tlst=list(np.arange(0,hh.shape[0],24))

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

#uavg=np.array(avgfile.u.isel(s_rho=lst).isel(ocean_time=tlst).values)

##########h=avgfile.groupby('ocean_time.month').mean().isel(s_rho=lst).u.values

#uavg = 0.5*(uavg[:,:,:,1:]+uavg[:,:,:,:-1])
lonu = 0.5*(lonu[:,1:]+lonu[:,:-1])
latu = 0.5*(latu[:,1:]+latu[:,:-1])


#uavg=uavg[:,:,1:-1,:]
lonu=lonu[1:-1,:]
latu=latu[1:-1,:]

##np.array(avgfile.u.isel(s_rho=lst).values)

#vavg=np.array(avgfile.v.isel(s_rho=lst).isel(ocean_time=tlst).values)
vavg=np.array(avgfile.groupby('ocean_time.month').mean().isel(s_rho=lst).v.values)


vavg = 0.5*(vavg[:,:,1:,:]+vavg[:,:,:-1,:])
latv = 0.5*(latv[1:,:]+latv[:-1,:])
lonv = 0.5*(lonv[1:,:]+lonv[:-1,:])

vavg=vavg[:,:,:,1:-1]
lonv=lonv[:,1:-1]
latv=latv[:,1:-1]


#tempavg=np.array(avgfile.temp.isel(s_rho=lst).values)
#tempavg=tempavg[:,:,1:-1,1:-1]


intu=np.zeros([vavg.shape[0],len(zc),x_roms.shape[0], x_roms.shape[1]])

intv=np.zeros([vavg.shape[0],len(zc),x_roms.shape[0], x_roms.shape[1]])

itemp=np.zeros([vavg.shape[0],len(zc),x_roms.shape[0], x_roms.shape[1]])

#zeta=avgfile['zeta'][:]

for j in range(intu.shape[2]):
  for k in range(intu.shape[3]):
    if (zr[-1,j,k] > zc.min()):
      zr[-1,j,k] = zc.min()

UNDEF=np.nan

for i in range(intu.shape[0]):
  for j in range(intu.shape[2]):
    for k in range(intu.shape[3]):
#      intu[i,:,j,k] = np.interp(-zc, -zr[:,j,k], uavg[i,:,j,k], right=UNDEF, left=UNDEF)
      intv[i,:,j,k] = np.interp(-zc, -zr[:,j,k], vavg[i,:,j,k], right=UNDEF, left=UNDEF)
#      itemp[i,:,j,k] = np.interp(-zc, -zr[:,j,k], tempavg[i,:,j,k], right=UNDEF, left=UNDEF)





lon=x_roms.copy()
lat=y_roms.copy()

intu[intu>100]=np.nan
u=np.ma.masked_invalid(intu)

intv[intv>100]=np.nan
v=np.ma.masked_invalid(intv)

itemp[itemp>100]=np.nan
temp=np.ma.masked_invalid(itemp)



lon  = x_roms.copy()
lat  = y_roms.copy()
#u2    = intu[2,:].copy()
v2    = intv[2,:].copy()


#v2[v2>0] = np.nan
#v2[v2<0] = np.nan

v2=np.ma.masked_invalid(v2)

vald=np.where(zc<3000)  
#vald=np.where(zc<400)  
vald=np.where((zc>400)&(zc<1000))  
row,col=np.where((lat<-20.5) & (lat>-23.5) & (lon<-39.4) & (lon>-40))
#row,col=np.where((lat<-20.5) & (lat>-23.5) & (lon<-39) & (lon>-40.74))  #23
row=np.unique(row)
col=np.unique(col)
col,row=np.meshgrid(col,row)


v2=np.squeeze(v2[vald,:])

v2=v2[:,row,col]

zi=zc[vald]

lon=lon[row,col]
lat=lat[row,col]

idxlat=np.argmin(np.abs(((np.abs(lat[:,0]) - 22))))


G,H=np.meshgrid(lon[idxlat,:], zi)
plt.pcolor(G,-H,v2[:,idxlat,:]);plt.colorbar();plt.show()




dx = np.diff(lon[idxlat,:], axis=0) * 111 * 1000# valid only for low latitudes!!!
aux= np.zeros([1]) 
aux[0]=dx[-1]
dx = np.concatenate( (dx, aux), axis=0) 
dz=abs(np.diff(zi, axis=0, append=0))
dx,dz=np.meshgrid(dx, dz)
transp = np.abs(dx) * np.abs(dz) * v2[:,idxlat,:]; transp = transp.sum()
print(transp / 1e6)

	

###############year

transpy=np.zeros(vavg.shape[0])


for i in range(len(transpy)):
  lon  = x_roms.copy()
  lat  = y_roms.copy()
  v2   = intv[i,:].copy()
  v2[v2>0] = np.nan
  #v2[v2<0] = np.nan
  v2=np.ma.masked_invalid(v2)
  vald=np.where(zc<400)  
  #vald=np.where((zc>400)&(zc<1400))  
  row,col=np.where((lat<-20.5) & (lat>-23.5) & (lon<-39.5) & (lon>-41.6))
  #row,col=np.where((lat<-20.5) & (lat>-23.5) & (lon<-39) & (lon>-40.74))  #23
  row=np.unique(row)
  col=np.unique(col)
  col,row=np.meshgrid(col,row)
  v2=np.squeeze(v2[vald,:])
  v2=v2[:,row,col]
  zi=zc[vald]
  lon=lon[row,col]
  lat=lat[row,col]
  idxlat=np.argmin(np.abs(((np.abs(lat[:,0]) - 22))))
  #G,H=np.meshgrid(lon[idxlat,:], zi)
  #plt.pcolor(G,-H,v2[:,idxlat,:]);plt.colorbar();plt.show()
  dx = np.diff(lon[idxlat,:], axis=0) * 111 * 1000# valid only for low latitudes!!!
  aux= np.zeros([1]) 
  aux[0]=dx[-1]
  dx = np.concatenate( (dx, aux), axis=0) 
  dz=abs(np.diff(zi, axis=0, append=0))
  dx,dz=np.meshgrid(dx, dz)
  transp = np.abs(dx) * np.abs(dz) * v2[:,idxlat,:]; transp = transp.sum()
  transpy[i] = transp/1e6

	
plt.pcolor(x_roms, y_roms, vavg[-1,-1,::]);plt.show()

[fig, ax] = plt.subplots(ncols=1, nrows=1, figsize=[8,3])
ax.plot(np.arange(1,13),transpy)
ax.fill_between(np.arange(1,13), transpy+transpy.std(), transpy-transpy.std(), alpha=0.5)
ax.set_xticks(np.arange(1,13))
ax.set_xticklabels(['jan', 'fev', 'mar','abr','mai','jun','jul','ago','set','out','nov','dez'])
ax.grid()
ax.set_ylabel('Volume Transport [Sv]')
plt.show()
#plt.savefig(path_fig+'temp_estatistica_mensal.png', dpi=200, bbox_inches='tight', transparent=False)



	

###############year

transpy=np.zeros(vavg.shape[0])


for i in range(len(transpy)):
  lon  = x_roms.copy()
  lat  = y_roms.copy()
  v2   = intv[i,:].copy()
  #v2[v2>0] = np.nan
  v2[v2<0] = np.nan
  v2=np.ma.masked_invalid(v2)
  #vald=np.where(zc<400)  
  vald=np.where((zc>500)&(zc<1500))  
  row,col=np.where((lat<-20.5) & (lat>-23.5) & (lon<-39.4) & (lon>-40))
  #row,col=np.where((lat<-20.5) & (lat>-23.5) & (lon<-39) & (lon>-40.74))  #23
  row=np.unique(row)
  col=np.unique(col)
  col,row=np.meshgrid(col,row)
  v2=np.squeeze(v2[vald,:])
  v2=v2[:,row,col]
  zi=zc[vald]
  lon=lon[row,col]
  lat=lat[row,col]
  idxlat=np.argmin(np.abs(((np.abs(lat[:,0]) - 22))))
  #G,H=np.meshgrid(lon[idxlat,:], zi)
  #plt.pcolor(G,-H,v2[:,idxlat,:]);plt.colorbar();plt.show()
  dx = np.diff(lon[idxlat,:], axis=0) * 111 * 1000# valid only for low latitudes!!!
  aux= np.zeros([1]) 
  aux[0]=dx[-1]
  dx = np.concatenate( (dx, aux), axis=0) 
  dz=abs(np.diff(zi, axis=0, append=0))
  dx,dz=np.meshgrid(dx, dz)
  transp = np.abs(dx) * np.abs(dz) * v2[:,idxlat,:]; transp = transp.sum()
  transpy[i] = transp/1e6


[fig, ax] = plt.subplots(ncols=1, nrows=1, figsize=[8,3])
ax.plot(np.arange(1,13),transpy)
ax.fill_between(np.arange(1,13), transpy+transpy.std(), transpy-transpy.std(), alpha=0.5)
ax.set_xticks(np.arange(1,13))
ax.set_xticklabels(['jan', 'fev', 'mar','abr','mai','jun','jul','ago','set','out','nov','dez'])
ax.grid()
ax.set_ylabel('Volume Transport [Sv]')
plt.show()
#plt.savefig(path_fig+'temp_estatistica_mensal.png', dpi=200, bbox_inches='tight', transparent=False)



#plot section

s_r=avgfile.s_rho.values
s_r[-1]=0

idxs=np.argmin(np.abs(((np.abs(y_roms[:,0]) - 20.8))))
#idxs=np.argmin(np.abs(((np.abs(y_roms[:,0]) - 22.7))))


#depm=np.zeros([len(s_r), len(x_roms[0,:])])
#for i in range(depm.shape[1]):
#  depm[:,i]=s_r*h_roms[idxs,i]

depm=-zr[:,idxs,:]

uu,oo=np.meshgrid(x_roms[idxs,:], np.arange(0,len(s_r)))

lind=depm[0,:]

levels = np.arange(-0.2,0.3,0.1)

vm=vavg[8:10,::]
vm=vm.mean(axis=0)

#vm=vavg[10,::]

vm=vavg.mean(axis=0)

[fig, ax] = plt.subplots(ncols=1, nrows=1, figsize=[7,5])
#yy=ax.pcolor(uu,depm, np.ma.masked_invalid(np.squeeze(vavg[7,:,idxs,:])), cmap=plt.get_cmap('seismic'), vmin=-0.4, vmax=0.4)
yy=ax.pcolormesh(uu,depm, np.ma.masked_invalid(np.squeeze(vm[:,idxs,:])), cmap=plt.get_cmap('seismic'), vmin=-0.4, vmax=0.4, shading='gouraud')
ax.fill_between(uu[0,:], lind, lind*0-5000, color='lightgrey', alpha=1)
cnt=ax.contour(uu,depm, np.ma.masked_invalid(np.squeeze(vm[:,idxs,:])),levels=levels,colors=('k'), linewidths=0.5)
plt.rcParams['contour.negative_linestyle'] = 'dashed'
ax.plot(uu[0,:], lind, color='black', linewidth=1.5)
#clabels=ax.clabel(cnt, inline=True, colors='k', fmt = '%1.1f',fontsize=6)
#[txt.set_path_effects([PathEffects.withStroke(linewidth=3, foreground='w')]) for txt in clabels]
cax = fig.add_axes([0.91, 0.3, 0.02, 0.38])
cbar=fig.colorbar(yy, shrink=0.8, extend='both', ax=ax, cax=cax)
cbar.ax.set_ylabel('Velocidade (m/s)', rotation=270)
cbar.ax.get_yaxis().labelpad = 11
cbar.ax.tick_params(labelsize=9)
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=12)
text.set_font_properties(font)
#ax.set_xlim([-41.5,-39.0])
ax.set_xlim([-40.5,-37.0])
ax.set_ylim([-2000,0])
ax.set_xlabel('Longitude (°)')
ax.set_ylabel('Profundidade (m)')
plt.savefig('21.png', dpi=200, bbox_inches='tight', transparent=False)
plt.show()
