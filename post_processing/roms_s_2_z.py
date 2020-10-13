from vgrid import s_coordinate, s_coordinate_2, s_coordinate_4
from netCDF4 import Dataset 
from scipy.interpolate import griddata, interp1d
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.pylab import *
from matplotlib import dates
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patheffects as PathEffects
import cmocean

   
     
avgfile=Dataset('HIS_FILE_20201012_5D0-20201019_5D0_fore_NEST_1_NEST_2.nc')

fname_grd = 'abc2.nc'

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
 
zc=np.array([4.940250e-01, 1.541375e+00, 2.645669e+00, 3.819495e+00,                                                                 5.078224e+00, 6.440614e+00, 7.929560e+00, 9.572997e+00,                                                                 1.140500e+01, 1.346714e+01, 1.581007e+01, 1.849556e+01,                                                                 2.159882e+01, 2.521141e+01, 2.944473e+01, 3.443415e+01,                                                                 4.034405e+01, 4.737369e+01, 5.576429e+01, 6.580727e+01,                                                                 7.785385e+01, 9.232607e+01, 1.097293e+02, 1.306660e+02,                                                                 1.558507e+02, 1.861256e+02, 2.224752e+02, 2.660403e+02,                                                                 3.181274e+02, 3.802130e+02, 4.539377e+02, 5.410889e+02,                                                                 6.435668e+02, 7.633331e+02, 9.023393e+02, 1.062440e+03,                                                                 1.245291e+03, 1.452251e+03, 1.684284e+03, 1.941893e+03,                                                                 2.225078e+03, 2.533336e+03, 2.865703e+03, 3.220820e+03,                                                                 3.597032e+03, 3.992484e+03, 4.405224e+03, 4.833291e+03,                                                               5.274784e+03, 5.727917e+03])
zc=zc[::-1]

uavg=avgfile['u_eastward'][:]

vavg=avgfile['v_northward'][:]

time=avgfile['ocean_time'][:]

tempavg=avgfile['temp'][:]

uavg=uavg[:,:]

vavg=vavg[:,:]

tempavg=tempavg[:,:]

time=time[:]/(24*60*60)

intu=np.zeros([uavg.shape[0],len(zc),x_roms.shape[0], x_roms.shape[1]])

intv=np.zeros([uavg.shape[0],len(zc),x_roms.shape[0], x_roms.shape[1]])

itemp=np.zeros([uavg.shape[0],len(zc),x_roms.shape[0], x_roms.shape[1]])

zeta=avgfile['zeta'][:]

for j in range(intu.shape[2]):
  for k in range(intu.shape[3]):
    if (zr[-1,j,k] > zc.min()):
      zr[-1,j,k] = zc.min()

UNDEF=np.nan

for i in range(intu.shape[0]):
  for j in range(intu.shape[2]):
    for k in range(intu.shape[3]):
#      I = interp1d(zr[:,j,k], uavg[i,:,j,k], kind='linear', bounds_error=False, assume_sorted=False)
#      intu[i,:,j,k] = I(zc)
#      I = interp1d(zr[:,j,k], vavg[i,:,j,k], kind='linear', bounds_error=False, assume_sorted=False)
#      intv[i,:,j,k] = I(zc)
#      I = interp1d(zr[:,j,k], tempavg[i,:,j,k], kind='linear', bounds_error=False, assume_sorted=False)
#      itemp[i,:,j,k] = I(zc)
      intu[i,:,j,k] = np.interp(-zc, -zr[:,j,k], uavg[i,:,j,k], right=UNDEF, left=UNDEF)
      intv[i,:,j,k] = np.interp(-zc, -zr[:,j,k], vavg[i,:,j,k], right=UNDEF, left=UNDEF)
      itemp[i,:,j,k] = np.interp(-zc, -zr[:,j,k], tempavg[i,:,j,k], right=UNDEF, left=UNDEF)

## test numpy
#UNDEF=np.nan
#      intu[i,:,j,k] = np.interp(-zc, -zr[:,j,k], uavg[i,:,j,k], right=UNDEF, left=UNDEF)
#      intv[i,:,j,k] = np.interp(-zc, -zr[:,j,k], vavg[i,:,j,k], right=UNDEF, left=UNDEF)
#      itemp[i,:,j,k] = np.interp(-zc, -zr[:,j,k], tempavg[i,:,j,k], right=UNDEF, left=UNDEF)
#


lon=x_roms.copy()
lat=y_roms.copy()

intu[intu>100]=np.nan
u=np.ma.masked_invalid(intu)

intv[intv>100]=np.nan
v=np.ma.masked_invalid(intv)

itemp[itemp>100]=np.nan
temp=np.ma.masked_invalid(itemp)


############################ u_ v eastward ja estão no TRUE East, positive) and northward (TRUE North, positive) d
#angle=grd['angle'][:]

#angle=np.tile(angle, (u.shape[0], u.shape[1],1,1))

#u=np.cos(angle)*u - np.sin(angle)*v

#v=np.cos(angle)*v + np.sin(angle)*u
#################################################


begindate=avgfile['ocean_time'].units[14:]
begindate=dates.datestr2num(begindate)
figdates=dates.num2date(begindate+time)
#figdates[tim].strftime("%A/%b - %-I %p" ) 
#figdates[tim].strftime("%m%d - %-I %p")

tim=-1
lay=-1

sp=2

up=u[tim,lay,:,:]

vp=v[tim,lay,:,:]

uu=up

vv=vp

val=np.sqrt((uu**2)+(vv**2))

a=plt.pcolor(x_roms,y_roms, np.squeeze(val),vmin=val.min(),vmax=1)
#a=plt.pcolor(x_roms,y_roms, np.squeeze(val),vmin=val.min(),vmax=val.max())
cnt=plt.contour(x_roms,y_roms,h_roms,levels=[200],colors=('red'))
plt.quiver(x_roms[0:-1:sp, 0:-1:sp],y_roms[0:-1:sp, 0:-1:sp],uu[0:-1:sp, 0:-1:sp], vv[0:-1:sp, 0:-1:sp])
#plt.text(-49,-28,figdates[tim].strftime("%m%d - %-I %p"), fontsize=12, fontweight='bold',
#        bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 5})
plt.text(-43,-25,figdates[tim].strftime("%m%d - %-I %p"), fontsize=12, fontweight='bold',
        bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 5})
plt.colorbar(a)
plt.show()



####################temp

############

begindate=avgfile['ocean_time'].units[14:]
begindate=dates.datestr2num(begindate)
figdates=dates.num2date(begindate+time)
#figdates[tim].strftime("%A/%b - %-I %p" ) 
#figdates[tim].strftime("%m%d - %-I %p")

tim=-1

lay=-1

a=temp[tim,lay,:]
aa=plt.pcolor(x_roms,y_roms, np.squeeze(a),vmin=21,vmax=27)
#aa=plt.pcolor(x_roms,y_roms, np.squeeze(a))
plt.colorbar(aa)
plt.text(-40,-18,figdates[tim].strftime("%m%d - %-I %p"), fontsize=12, fontweight='bold',
        bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 5})
plt.show()



#plt.pcolor(x_roms,y_roms, h_roms,vmax=200);plt.colorbar();plt.show()

################WIND



angle=grd['angle'][:]
uwind=avgfile['Uwind'][:]
vwind=avgfile['Vwind'][:]



sp=3
tim=1

uu = uwind[tim,:]
vv = vwind[tim,:]

uout=np.cos(angle)*uu - np.sin(angle)*vv

vout=np.cos(angle)*vv + np.sin(angle)*uu

#uout=uu

#vout=vv

val=np.sqrt((uout[:,:]**2)+(vout[:,:]**2))

a=plt.pcolor(x_roms,y_roms, np.squeeze(val),vmin=0, vmax=13)
#a=plt.pcolor(x_roms,y_roms, np.squeeze(val),vmin=0)
plt.colorbar(a)
plt.quiver(x_roms[0:-1:sp, 0:-1:sp],y_roms[0:-1:sp, 0:-1:sp],uout[0:-1:sp, 0:-1:sp], vout[0:-1:sp, 0:-1:sp])
plt.text(-40,-18,figdates[tim].strftime("%m%d - %-I %p"), fontsize=12, fontweight='bold',
        bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 5})
plt.show()




tim=0

plt.pcolor(lon,lat,zeta[tim,:]);plt.colorbar();plt.show()





tim=10

lay=-1

sp=2

val=np.sqrt((u[tim,lay,:,:]**2)+(v[tim,lay,:,:]**2))

a=plt.pcolor(x_roms,y_roms, np.squeeze(val),vmin=val.min(),vmax=val.max())
#a=plt.pcolor(x_roms,y_roms, np.squeeze(val),vmin=val.min(),vmax=val.max())
plt.quiver(x_roms[0:-1:sp, 0:-1:sp],y_roms[0:-1:sp, 0:-1:sp],u[tim,lay,0:-1:sp, 0:-1:sp], v[tim,lay,0:-1:sp, 0:-1:sp])
plt.colorbar(a)
plt.show()



#########################################make gif
begindate=avgfile['ocean_time'].units[14:]
begindate=dates.datestr2num(begindate)
figdates=dates.num2date(begindate+time)
#figdates[tim].strftime("%A/%b - %-I %p" ) 
#figdates[tim].strftime("%m%d - %-I %p")

sp=4
lay=-1

rj=(-42.7,-22.5)
mg=(-43.5,-19.5)
es=(-40.5,-19.5)
datest=(-50,-22.5)

#nlim=-18
#slim=-27
#wlim=-49
#elim=-37

nlim=-20.7
slim=lat.min()
wlim=lon.min()
elim=lon.max()

pathfig='/home/fernando/roms/src/Projects/operational/figure/SUL/'

for tim in range(u.shape[0]):
#for tim in [0]:
  val=np.sqrt((u[tim,lay,:,:]**2)+(v[tim,lay,:,:]**2))
  fig, axs = plt.subplots(figsize=(15,15))
  map = Basemap(projection='merc', llcrnrlat=slim, urcrnrlat=nlim,llcrnrlon=wlim, urcrnrlon=elim,resolution='l')
  map.drawcoastlines(linewidth=0.25)
  map.drawcountries(linewidth=0.25)
  map.fillcontinents(color='#ddaa66',lake_color='aqua')
  parallels = np.arange(-40,-15,2)
  map.drawparallels(parallels,labels=[1,0,0,1], linewidth=0.0)
  meridians = np.arange(-55,-30,2)
  map.drawmeridians(meridians,labels=[1,0,0,1], linewidth=0.0)
  map.readshapefile('/mnt/c/Users/Fernando/Desktop/shape_unid_fed/lim_unidade_federacao_a', 'lim_unidade_federacao_a')
  x, y = map(lon,lat)
  a=axs.pcolor(x,y, np.squeeze(val),vmin=0,vmax=0.9)
  axs.quiver(x[0:-1:sp, 0:-1:sp],y[0:-1:sp, 0:-1:sp],u[tim,lay,0:-1:sp, 0:-1:sp], v[tim,lay,0:-1:sp, 0:-1:sp])
#  axs.text(map(rj[0], rj[1])[0],map(rj[0], rj[1])[1],'RJ', fontsize=12, fontweight='bold',
#        bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 4})
#  axs.text(map(mg[0], mg[1])[0],map(mg[0], mg[1])[1],'MG', fontsize=14, fontweight='bold',
#        bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 5})
#  axs.text(map(es[0], es[1])[0],map(es[0], es[1])[1],'ES', fontsize=14, fontweight='bold',
#        bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 5})
  axs.text(map(datest[0], datest[1])[0],map(datest[0], datest[1])[1],figdates[tim].strftime("%m%d - %-I %p"), fontsize=18, fontweight='bold',
        bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 5})
  divider = make_axes_locatable(axs)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  cbar=fig.colorbar(a, shrink=0.8, extend='both', ax=axs, cax=cax)
  cbar.ax.set_ylabel('Velocity(m/s)', rotation=270)
  cbar.ax.get_yaxis().labelpad = 20
  text = cbar.ax.yaxis.label
  font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=15)
  text.set_font_properties(font)
  plt.savefig(pathfig+str(tim)+'father.png')

convert -resize 768x576 -delay 35 -loop 0 `ls -v *.png` myimage.gif   #maior delay mais lento

ffmpeg -i myimage.gif -movflags faststart -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" forecast.mp4

cp forecast.mp4 /mnt/c/Users/Fernando/Desktop/ 

mv myimage.gif /mnt/c/Users/Fernando/Desktop/ 



convert -resize 1000x800 -delay 35 -loop 0 `ls -v *.png` myimage2.gif   #maior delay mais lento

ffmpeg -i myimage2.gif -movflags faststart -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" forecast.mp4

cp forecast.mp4 /mnt/c/Users/Fernando/Desktop




##zoom

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

tim=-1
lay=-1

sp=5

up=u[tim,lay,:,:]

vp=v[tim,lay,:,:]

uu=up

vv=vp

val=np.sqrt((uu**2)+(vv**2))

fig, ax1 = plt.subplots()

a=ax1.pcolor(x_roms,y_roms, np.squeeze(val),vmin=val.min(),vmax=0.80)
#a=plt.pcolor(x_roms,y_roms, np.squeeze(val),vmin=val.min(),vmax=val.max())
cnt=ax1.contour(x_roms,y_roms,h_roms,levels=[200],colors=('red'))
ax1.quiver(x_roms[0:-1:sp, 0:-1:sp],y_roms[0:-1:sp, 0:-1:sp],uu[0:-1:sp, 0:-1:sp], vv[0:-1:sp, 0:-1:sp])
#plt.text(-49,-28,figdates[tim].strftime("%m%d - %-I %p"), fontsize=12, fontweight='bold',
#        bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 5})
plt.text(-44,-21,figdates[tim].strftime("%m%d - %-I %p"), fontsize=12, fontweight='bold',
        bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 5})
plt.colorbar(a)
#plt.show()



axins = zoomed_inset_axes(ax1, 2, loc=1) # zoom-factor: 2.5, location: upper-left
axins.plot(x_roms,y_roms)
x1, x2, y1, y2 = -44, -43, -24, -23.5 # specify the limits
axins.set_xlim(x1, x2) # apply the x-limits
axins.set_ylim(y1, y2) # apply the y-limits
axins.tick_params(axis="y",direction="in", pad=-22)
axins.tick_params(axis="x",direction="in", pad=-15)
plt.yticks(visible=False)
plt.xticks(visible=False)
mark_inset(ax1, axins, loc1=2, loc2=1, fc="none", ec="0.5")
plt.show()



################################## zooooooooooooommmm

from vgrid import s_coordinate, s_coordinate_2, s_coordinate_4
from netCDF4 import Dataset 
from scipy.interpolate import griddata, interp1d
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.pylab import *
from matplotlib import dates
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patheffects as PathEffects
import cmocean
     
avgfile=Dataset('HIS_FILE_20201012_5D0-20201019_5D0_fore.nc')

fname_grd = 'azul_grd2.nc'

## Load ROMS grid.
grd = Dataset(fname_grd)
x_roms = grd.variables['lon_rho'][:]
y_roms = grd.variables['lat_rho'][:]
msk_roms = grd.variables['mask_rho'][:]
msk_romsv = grd.variables['mask_v'][:]
h_roms = grd.variables['h'][:]
h_roms1=h_roms.copy()
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
 
zc=np.array([4.940250e-01, 1.541375e+00, 2.645669e+00, 3.819495e+00,                                                                 5.078224e+00, 6.440614e+00, 7.929560e+00, 9.572997e+00,                                                                 1.140500e+01, 1.346714e+01, 1.581007e+01, 1.849556e+01,                                                                 2.159882e+01, 2.521141e+01, 2.944473e+01, 3.443415e+01,                                                                 4.034405e+01, 4.737369e+01, 5.576429e+01, 6.580727e+01,                                                                 7.785385e+01, 9.232607e+01, 1.097293e+02, 1.306660e+02,                                                                 1.558507e+02, 1.861256e+02, 2.224752e+02, 2.660403e+02,                                                                 3.181274e+02, 3.802130e+02, 4.539377e+02, 5.410889e+02,                                                                 6.435668e+02, 7.633331e+02, 9.023393e+02, 1.062440e+03,                                                                 1.245291e+03, 1.452251e+03, 1.684284e+03, 1.941893e+03,                                                                 2.225078e+03, 2.533336e+03, 2.865703e+03, 3.220820e+03,                                                                 3.597032e+03, 3.992484e+03, 4.405224e+03, 4.833291e+03,                                                               5.274784e+03, 5.727917e+03])
zc=zc[::-1]

uavg=avgfile['u_eastward'][:]
vavg=avgfile['v_northward'][:]
time=avgfile['ocean_time'][:]
tempavg=avgfile['temp'][:]
uavg=uavg[:,:]
vavg=vavg[:,:]
tempavg=tempavg[:,:]
time=time[:]/(24*60*60)
intu=np.zeros([uavg.shape[0],len(zc),x_roms.shape[0], x_roms.shape[1]])
intv=np.zeros([uavg.shape[0],len(zc),x_roms.shape[0], x_roms.shape[1]])
itemp=np.zeros([uavg.shape[0],len(zc),x_roms.shape[0], x_roms.shape[1]])
zeta=avgfile['zeta'][:]
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
#      itemp[i,:,j,k] = np.interp(-zc, -zr[:,j,k], tempavg[i,:,j,k], right=UNDEF, left=UNDEF)



lon1=x_roms.copy()
lat1=y_roms.copy()
intu[intu>100]=np.nan
u1=np.ma.masked_invalid(intu)
intv[intv>100]=np.nan
v1=np.ma.masked_invalid(intv)
#itemp[itemp>100]=np.nan
#temp=np.ma.masked_invalid(itemp)



     
avgfile=Dataset('HIS_FILE_20201012_5D0-20201019_5D0_fore_NEST_1.nc')

fname_grd = 'abc1.nc'

## Load ROMS grid.
grd = Dataset(fname_grd)
x_roms = grd.variables['lon_rho'][:]
y_roms = grd.variables['lat_rho'][:]
msk_roms = grd.variables['mask_rho'][:]
msk_romsv = grd.variables['mask_v'][:]
h_roms = grd.variables['h'][:]
h_roms2=h_roms.copy()
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
 
zc=np.array([4.940250e-01, 1.541375e+00, 2.645669e+00, 3.819495e+00,                                                                 5.078224e+00, 6.440614e+00, 7.929560e+00, 9.572997e+00,                                                                 1.140500e+01, 1.346714e+01, 1.581007e+01, 1.849556e+01,                                                                 2.159882e+01, 2.521141e+01, 2.944473e+01, 3.443415e+01,                                                                 4.034405e+01, 4.737369e+01, 5.576429e+01, 6.580727e+01,                                                                 7.785385e+01, 9.232607e+01, 1.097293e+02, 1.306660e+02,                                                                 1.558507e+02, 1.861256e+02, 2.224752e+02, 2.660403e+02,                                                                 3.181274e+02, 3.802130e+02, 4.539377e+02, 5.410889e+02,                                                                 6.435668e+02, 7.633331e+02, 9.023393e+02, 1.062440e+03,                                                                 1.245291e+03, 1.452251e+03, 1.684284e+03, 1.941893e+03,                                                                 2.225078e+03, 2.533336e+03, 2.865703e+03, 3.220820e+03,                                                                 3.597032e+03, 3.992484e+03, 4.405224e+03, 4.833291e+03,                                                               5.274784e+03, 5.727917e+03])
zc=zc[::-1]

uavg=avgfile['u_eastward'][:]
vavg=avgfile['v_northward'][:]
time=avgfile['ocean_time'][:]
tempavg=avgfile['temp'][:]
uavg=uavg[:,:]
vavg=vavg[:,:]
tempavg=tempavg[:,:]
time=time[:]/(24*60*60)
intu=np.zeros([uavg.shape[0],len(zc),x_roms.shape[0], x_roms.shape[1]])
intv=np.zeros([uavg.shape[0],len(zc),x_roms.shape[0], x_roms.shape[1]])
itemp=np.zeros([uavg.shape[0],len(zc),x_roms.shape[0], x_roms.shape[1]])
zeta=avgfile['zeta'][:]
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
#      itemp[i,:,j,k] = np.interp(-zc, -zr[:,j,k], tempavg[i,:,j,k], right=UNDEF, left=UNDEF)



lon2=x_roms.copy()
lat2=y_roms.copy()
intu[intu>100]=np.nan
u2=np.ma.masked_invalid(intu)
intv[intv>100]=np.nan
v2=np.ma.masked_invalid(intv)
#itemp[itemp>100]=np.nan
#temp2=np.ma.masked_invalid(itemp)





     
avgfile=Dataset('HIS_FILE_20201012_5D0-20201019_5D0_fore_NEST_1_NEST_2.nc')

fname_grd = 'abc2.nc'

## Load ROMS grid.
grd = Dataset(fname_grd)
x_roms = grd.variables['lon_rho'][:]
y_roms = grd.variables['lat_rho'][:]
msk_roms = grd.variables['mask_rho'][:]
msk_romsv = grd.variables['mask_v'][:]
h_roms = grd.variables['h'][:]
h_roms3=h_roms.copy()
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
 
zc=np.array([4.940250e-01, 1.541375e+00, 2.645669e+00, 3.819495e+00,                                                                 5.078224e+00, 6.440614e+00, 7.929560e+00, 9.572997e+00,                                                                 1.140500e+01, 1.346714e+01, 1.581007e+01, 1.849556e+01,                                                                 2.159882e+01, 2.521141e+01, 2.944473e+01, 3.443415e+01,                                                                 4.034405e+01, 4.737369e+01, 5.576429e+01, 6.580727e+01,                                                                 7.785385e+01, 9.232607e+01, 1.097293e+02, 1.306660e+02,                                                                 1.558507e+02, 1.861256e+02, 2.224752e+02, 2.660403e+02,                                                                 3.181274e+02, 3.802130e+02, 4.539377e+02, 5.410889e+02,                                                                 6.435668e+02, 7.633331e+02, 9.023393e+02, 1.062440e+03,                                                                 1.245291e+03, 1.452251e+03, 1.684284e+03, 1.941893e+03,                                                                 2.225078e+03, 2.533336e+03, 2.865703e+03, 3.220820e+03,                                                                 3.597032e+03, 3.992484e+03, 4.405224e+03, 4.833291e+03,                                                               5.274784e+03, 5.727917e+03])
zc=zc[::-1]

uavg=avgfile['u_eastward'][:]
vavg=avgfile['v_northward'][:]
time=avgfile['ocean_time'][:]
tempavg=avgfile['temp'][:]
uavg=uavg[:,:]
vavg=vavg[:,:]
tempavg=tempavg[:,:]
time=time[:]/(24*60*60)
intu=np.zeros([uavg.shape[0],len(zc),x_roms.shape[0], x_roms.shape[1]])
intv=np.zeros([uavg.shape[0],len(zc),x_roms.shape[0], x_roms.shape[1]])
itemp=np.zeros([uavg.shape[0],len(zc),x_roms.shape[0], x_roms.shape[1]])
zeta=avgfile['zeta'][:]
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
#      itemp[i,:,j,k] = np.interp(-zc, -zr[:,j,k], tempavg[i,:,j,k], right=UNDEF, left=UNDEF)



lon3=x_roms.copy()
lat3=y_roms.copy()
intu[intu>100]=np.nan
u3=np.ma.masked_invalid(intu)
intv[intv>100]=np.nan
v3=np.ma.masked_invalid(intv)
#itemp[itemp>100]=np.nan
#temp2=np.ma.masked_invalid(itemp)



from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

tim=-1
lay=-1

pathfig='/home/fernando/roms/src/Projects/operational/figure/proo/'
for tim in range(u1.shape[0]):
#for tim in [0]:
  sp=4
  up=u1[tim,lay,:,:]
  vp=v1[tim,lay,:,:]
  uu=up
  vv=vp
  val=np.sqrt((uu**2)+(vv**2))
  begindate=avgfile['ocean_time'].units[14:]
  begindate=dates.datestr2num(begindate)
  figdates=dates.num2date(begindate+time)
  nlim=lat1.max()
  slim=lat1.min()
  wlim=lon1.min()
  elim=lon1.max()
#
  fig, ax1 = plt.subplots(figsize=(8,8))
  map = Basemap(projection='cyl', llcrnrlat=slim, urcrnrlat=nlim,llcrnrlon=wlim, urcrnrlon=elim,resolution='l')
  map.drawcoastlines(linewidth=0.25)
  map.drawcountries(linewidth=0.25)
  map.fillcontinents(color='silver',lake_color='white')
  parallels = np.arange(-40,-15,3)
  map.drawparallels(parallels,labels=[1,0,0,1], linewidth=0.0)
  meridians = np.arange(-55,-15,3)
  map.drawmeridians(meridians,labels=[1,0,0,1], linewidth=0.0)
  map.readshapefile('/mnt/c/Users/Fernando/Desktop/shape_unid_fed/lim_unidade_federacao_a', 'lim_unidade_federacao_a')
#  map.drawmapscale(-39., -24.5, -39, -24.5, 100, barstyle='fancy', fontsize = 8, yoffset=5000)
  x_roms, y_roms = map(lon1,lat1)
#
  a=ax1.pcolor(x_roms,y_roms, np.squeeze(val),vmin=val.min(),vmax=1.0)
  ax1.quiver(x_roms[0:-1:sp, 0:-1:sp],y_roms[0:-1:sp, 0:-1:sp],uu[0:-1:sp, 0:-1:sp], vv[0:-1:sp, 0:-1:sp])
  ax1.contourf(x_roms,y_roms,h_roms1,levels=[0,20],colors=('gainsboro'))
  plt.text(-49.5,-17,figdates[tim].strftime("%b/%d - %-I %p"), fontsize=10, fontweight='bold',
        bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 5})
  divider = make_axes_locatable(ax1)
  #cax = divider.append_axes("right", size="5%", pad=0.05)
  cax = fig.add_axes([0.91, 0.3, 0.02, 0.38])
  cbar=fig.colorbar(a, shrink=0.8, extend='both', ax=ax1, cax=cax)
  cbar.ax.set_ylabel('Velocity(m/s)', rotation=270)
  cbar.ax.get_yaxis().labelpad = 11
  text = cbar.ax.yaxis.label
  font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=11)
  text.set_font_properties(font)
#
#
  sp=3
  up=u2[tim,lay,:,:]
  vp=v2[tim,lay,:,:]
  uu=up
  vv=vp
  val=np.sqrt((uu**2)+(vv**2))
#
#
  left, bottom, width, height = [0.50, 0.60, 0.35, 0.35]
  axins = fig.add_axes([left, bottom, width, height])
#
#
  nlim=lat2.max()
  slim=lat2.min()
  wlim=lon2.min()
  elim=lon2.max()
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
#  map.drawmapscale(-39., -24.5, -39, -24.5, 100, barstyle='fancy', fontsize = 8, yoffset=5000)
  x_roms2, y_roms2 = map(lon2,lat2)
#
#
  axins.pcolor(x_roms2,y_roms2, np.squeeze(val),vmin=val.min(),vmax=1.0)
  axins.contourf(x_roms2,y_roms2,h_roms2,levels=[0,20],colors=('gainsboro'))
  axins.quiver(x_roms2[0:-1:sp, 0:-1:sp],y_roms2[0:-1:sp, 0:-1:sp],uu[0:-1:sp, 0:-1:sp], vv[0:-1:sp, 0:-1:sp])
#
  mark_inset(ax1, axins, loc1=2, loc2=4, fc="none", ec="black")
#
#
  sp=3
  up=u3[tim,lay,:,:]
  vp=v3[tim,lay,:,:]
  uu=up
  vv=vp
  val=np.sqrt((uu**2)+(vv**2))
#
#
  left, bottom, width, height = [0.50, 0.05, 0.35, 0.35]
  axins2 = fig.add_axes([left, bottom, width, height])
#
#
  nlim=lat3.max()
  slim=lat3.min()
  wlim=lon3.min()
  elim=lon3.max()
#
  map = Basemap(projection='cyl', llcrnrlat=slim, urcrnrlat=nlim,llcrnrlon=wlim, urcrnrlon=elim,resolution='l')
  map.drawcoastlines(linewidth=0.25)
  map.drawcountries(linewidth=0.25)
  map.fillcontinents(color='silver',lake_color='white')
  parallels = np.arange(-40,-15,1)
  map.drawparallels(parallels,labels=[0,0,0,0], linewidth=0.0)
  meridians = np.arange(-55,-30,1)
  map.drawmeridians(meridians,labels=[0,0,0,0], linewidth=0.0)
  map.readshapefile('/mnt/c/Users/Fernando/Desktop/shape_unid_fed/lim_unidade_federacao_a', 'lim_unidade_federacao_a')
#map.drawmapscale(-39., -24.5, -39, -24.5, 100, barstyle='fancy', fontsize = 8, yoffset=5000)
  x_roms3, y_roms3 = map(lon3,lat3)
#
#
  axins2.pcolor(x_roms3,y_roms3, np.squeeze(val),vmin=val.min(),vmax=1.0)
  axins2.quiver(x_roms3[0:-1:sp, 0:-1:sp],y_roms3[0:-1:sp, 0:-1:sp],uu[0:-1:sp, 0:-1:sp], vv[0:-1:sp, 0:-1:sp])
  axins2.contourf(x_roms3,y_roms3,h_roms3,levels=[0,20],colors=('gainsboro'))
#
  mark_inset(ax1, axins2, loc1=3, loc2=1, fc="none", ec="black")
#
  for axis in ['top','bottom','left','right']:
    ax1.spines[axis].set_linewidth(1.5)
    ax1.spines[axis].set_color('black')
#
##plt.tight_layout()
  plt.savefig(pathfig+str(tim)+'ope.png')

convert -resize 1000x800 -delay 45 -loop 0 `ls -v *.png` myimage2.gif   #maior delay mais lento

ffmpeg -i myimage2.gif -movflags faststart -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" forecast_p.mp4

mv forecast_p.mp4 /mnt/c/Users/Fernando/Desktop
cp myimage2.gif /mnt/c/Users/Fernando/Desktop



#######################################################

##sao tome eddy

begindate=avgfile['ocean_time'].units[14:]
begindate=dates.datestr2num(begindate)
figdates=dates.num2date(begindate+time)
#figdates[tim].strftime("%A/%b - %-I %p" ) 
#figdates[tim].strftime("%m%d - %-I %p")

tim=-1
lay=-1

sp=3

up=u[tim,lay,:,:]

vp=v[tim,lay,:,:]

uu=up

vv=vp

val=np.sqrt((uu**2)+(vv**2))

nlim=lat.max()-0.5
slim=lat.min()
wlim=lon.min()
elim=lon.max()


fig, axs = plt.subplots(figsize=(10,10))
map = Basemap(projection='merc', llcrnrlat=slim, urcrnrlat=nlim,llcrnrlon=wlim, urcrnrlon=elim,resolution='l')
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.fillcontinents(color='tan',lake_color='white')
parallels = np.arange(-40,-15,1)
map.drawparallels(parallels,labels=[1,0,0,1], linewidth=0.0)
meridians = np.arange(-55,-30,1)
map.drawmeridians(meridians,labels=[1,0,0,1], linewidth=0.0)
map.readshapefile('/mnt/c/Users/Fernando/Desktop/shape_unid_fed/lim_unidade_federacao_a', 'lim_unidade_federacao_a')
map.drawmapscale(-39., -24.5, -39, -24.5, 100, barstyle='fancy', fontsize = 8, yoffset=5000)
x, y = map(lon,lat)


a=axs.pcolor(x,y, np.squeeze(val),vmin=val.min(),vmax=1.0)
#a=plt.pcolor(x_roms,y_roms, np.squeeze(val),vmin=val.min(),vmax=val.max())
cnt=axs.contour(x,y,h_roms,levels=[30,200],colors=('k'), linewidths=2.0)
axs.contourf(x,y,h_roms,levels=[0,30],colors=('lightyellow'))
#axs.contourf(x,y,h_roms,levels=[0,30],cmap=cmocean.cm.cmap_d['deep'])
clabels=axs.clabel(cnt, inline=False, colors='k', fmt = '%2.1d',fontsize=8)
[txt.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')]) for txt in clabels]
axs.quiver(x[0:-1:sp, 0:-1:sp],y[0:-1:sp, 0:-1:sp],uu[0:-1:sp, 0:-1:sp], vv[0:-1:sp, 0:-1:sp])
#plt.text(-49,-28,figdates[tim].strftime("%m%d - %-I %p"), fontsize=12, fontweight='bold',
#        bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 5})
#plt.text(-42,-22,figdates[tim].strftime("%m%d - %-I %p"), fontsize=12, fontweight='bold',
#        bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 5})
divider = make_axes_locatable(axs)
#cax = divider.append_axes("right", size="5%", pad=0.05)
cax = fig.add_axes([0.83, 0.3, 0.02, 0.38])
cbar=fig.colorbar(a, shrink=0.8, extend='both', ax=axs, cax=cax)
cbar.ax.set_ylabel('Velocity(m/s)', rotation=270)
cbar.ax.get_yaxis().labelpad = 20
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=13)
text.set_font_properties(font)


left, bottom, width, height = [0.20, 0.7, 0.15, 0.15]
ax3 = fig.add_axes([left, bottom, width, height])

map = Basemap(projection='ortho',lat_0=-20,lon_0=-50,resolution='l', ax=ax3)
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.fillcontinents(color='tan',lake_color='aqua')
map.drawmapboundary(fill_color='powderblue')
#map.drawmeridians(np.arange(0,360,30))
#map.drawparallels(np.arange(-90,90,30))
#x, y = map(lon,lat)
#ax2.plot(x,y, 'blue')


x1,y1 = map(-50,-30)
x2,y2 = map(-50,-15)
x3,y3 = map(-30,-15)
x4,y4 = map(-30,-30)
poly = Polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4)],facecolor='blue',edgecolor='black',linewidth=1, alpha=0.8)
ax3.add_patch(poly)

#plt.tight_layout()

for axis in ['top','bottom','left','right']:
  axs.spines[axis].set_linewidth(1.5)
  axs.spines[axis].set_color('black')
  

plt.show()

