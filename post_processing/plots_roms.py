import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset as dat
from vgrid import s_coordinate, s_coordinate_2, s_coordinate_4
from matplotlib.pylab import *



file=dat('20200421.5D0_NEST_1_bry.nc')
#file1=dat('prooceano_myocean_mais_pontos_12134_bry.nc')

temp=file['v_south'][::]

#temp1=file1['v_north'][::]

#plt.pcolor(temp[0,:,:]);plt.colorbar();plt.show()

#plt.pcolor(temp1[7,:,:]);plt.colorbar();plt.show()


theta_b = 0.4
theta_s = 5.0
tcline = 100.
klevels = 40
Vtransform = 2
Vstretching = 4
fname_grd='azul_son_case1_newBAT_2_Mcu.nc'
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

a,b=np.meshgrid(x_roms[0,:],zr[:,0,0])

fig, ax1 = plt.subplots()
AB=ax1.pcolor(a, zr[:,0,:], temp[-1,::], vmin=-0.8, vmax=0.8, cmap = 'seismic')
#AB=ax1.pcolor(a, zr[:,0,:], temp[0,::])
ax1.axis('tight')
cbar=fig.colorbar(AB, shrink=0.9, extend='both')
cbar.ax.set_ylabel('v(m/s)', rotation=270)
cbar.ax.get_yaxis().labelpad = 20
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=15)
text.set_font_properties(font)
#plt.savefig('bry_south.png')
plt.ylim(-20,3)


plt.show()

a,b=np.meshgrid(x_roms[0,:],zr[:,0,0])

fig, ax1 = plt.subplots()
AB=ax1.pcolor(a, zr[:,-1,:], temp1[0,::], vmin=-0.8, vmax=0.8, cmap = 'seismic')
#AB=ax1.pcolor(a, zr[:,0,:], temp[0,::])
ax1.axis('tight')
cbar=fig.colorbar(AB, shrink=0.9, extend='both')
cbar.ax.set_ylabel('v(m/s)', rotation=270)
cbar.ax.get_yaxis().labelpad = 20
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=15)
text.set_font_properties(font)
#plt.savefig('bry_south.png')
plt.ylim(-5000,0)


plt.show()

#DIFFERENCE NORTH/SOUTH

tempdif = temp - temp1

a,b=np.meshgrid(x_roms[0,:],zr[:,0,0])

fig, ax1 = plt.subplots()
AB=ax1.pcolor(a, zr[:,-1,:], tempdif[-1,::], vmin=-0.3, vmax=0.3, cmap = 'seismic')
#AB=ax1.pcolor(a, zr[:,0,:], temp[0,::])
ax1.axis('tight')
cbar=fig.colorbar(AB, shrink=0.9, extend='both')
cbar.ax.set_ylabel('v(m/s)', rotation=270)
cbar.ax.get_yaxis().labelpad = 20
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=15)
text.set_font_properties(font)
#plt.savefig('bry_south.png')
plt.ylim(-5000,0)


plt.show()



#EAST


a,b=np.meshgrid(y_roms[:,0],zr[:,0,0])

fig, ax1 = plt.subplots()
AB=ax1.pcolor(a, zr[:,:,-1], temp[1,::], cmap = 'RdBu')
#AB=ax1.pcolor(a, zr[:,:,-1], temp[0,::])
ax1.axis('tight')
cbar=fig.colorbar(AB, shrink=0.9, extend='both')
cbar.ax.set_ylabel('v($m^{2}$)', rotation=270)
cbar.ax.get_yaxis().labelpad = 20
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=15)
text.set_font_properties(font)

plt.show()



a,b=np.meshgrid(y_roms[:,0],zr[:,0,0])

fig, ax1 = plt.subplots()
AB=ax1.pcolor(a, zr[:,:,-1], temp1[-1,::], vmin=-0.15, vmax=0.15, cmap = 'RdBu')
#AB=ax1.pcolor(a, zr[:,:,-1], temp[0,::])
ax1.axis('tight')
cbar=fig.colorbar(AB, shrink=0.9, extend='both')
cbar.ax.set_ylabel('v($m^{2}$)', rotation=270)
cbar.ax.get_yaxis().labelpad = 20
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=15)
text.set_font_properties(font)

plt.show()



#DIFFERENCE EAST

tempdif = temp - temp1

a,b=np.meshgrid(y_roms[:,0],zr[:,0,0])

fig, ax1 = plt.subplots()
AB=ax1.pcolor(a, zr[:,:,-1], tempdif[-1,::], vmin=-0.5, vmax=0.5, cmap = 'seismic')
#AB=ax1.pcolor(a, zr[:,:,-1], temp[0,::])
ax1.axis('tight')
cbar=fig.colorbar(AB, shrink=0.9, extend='both')
cbar.ax.set_ylabel('v($m^{2}$)', rotation=270)
cbar.ax.get_yaxis().labelpad = 20
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=15)
text.set_font_properties(font)

plt.show()





############################# INITIAL FILE
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset as dat
from vgrid import s_coordinate, s_coordinate_2, s_coordinate_4
from matplotlib.pylab import *

file=dat('20200726_ini.nc')



theta_b = 0.4
theta_s = 5.0
tcline = 100.
klevels = 40
Vtransform = 2
Vstretching = 4

fname_grd = 'grid_rotated_SUL_2_smaler.nc'

grd = Dataset(fname_grd)
x_roms = grd.variables['lon_rho'][:]
y_roms = grd.variables['lat_rho'][:]
msk_roms = grd.variables['mask_rho'][:]
msk_romsv = grd.variables['mask_v'][:]
h_roms = grd.variables['h'][:]
sh2 = msk_roms.shape
etamax, ximax = sh2
angle=grd['angle'][:]

scoord = s_coordinate_4(h_roms, theta_b, theta_s, tcline, klevels)
zr = scoord.z_r[:]


tempini=file['u'][0,:,-1,  :]



a,b=np.meshgrid(x_roms[-1,:],zr[:,-1,0])

fig, ax1 = plt.subplots()
AB=ax1.pcolor(a, zr[:,-1,:], np.squeeze(tempini), vmin=tempini.min(), vmax=tempini.max(), cmap = 'seismic')
#AB=ax1.pcolor(a, zr[:,0,:], temp[0,::])
ax1.axis('tight')
cbar=fig.colorbar(AB, shrink=0.9, extend='both')
cbar.ax.set_ylabel('v(m/s)', rotation=270)
cbar.ax.get_yaxis().labelpad = 20
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=15)
text.set_font_properties(font)
#plt.savefig('bry_south.png')
plt.ylim(-5000,0)
plt.show()


angle=angle[:-1,:-1]



tim=-1
lay=-1

vp=file['v'][tim,lay,:,  :-1]

up=file['u'][tim,lay,:-1,  :]


u=np.cos(angle)*up - np.sin(angle)*vp

v=np.cos(angle)*vp + np.sin(angle)*up

val=np.sqrt((u**2)+(v**2))


plt.pcolor(x_roms, y_roms, val, vmin=val.min(),vmax=1.1);plt.colorbar();plt.show()


plt.pcolor(x_roms, y_roms, v);plt.colorbar();plt.show()



plt.pcolor(x_roms, y_roms, u,  vmax=0.5, vmin=-0.7);plt.colorbar();plt.show()


############################# CLIMATOLOGY FILE
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset as dat
from vgrid import s_coordinate, s_coordinate_2, s_coordinate_4
from matplotlib import *

file=Dataset('grid_rotated_clm_clm.nc')
tempini=file['v'][::]

time=file['temp_time'][:]/(24*60*60)
begindate=file['temp_time'].units[14:]
begindate=dates.datestr2num(begindate)
figdates=dates.num2date(begindate+time)


theta_b = 0.4
theta_s = 5.0
tcline = 100.
klevels = 40
Vtransform = 2
Vstretching = 4

fname_grd='grid_rotated_SUL_2_NEST_2_lp.nc'

grd = Dataset(fname_grd)
x_roms = grd.variables['lon_rho'][:]
y_roms = grd.variables['lat_rho'][:]
msk_roms = grd.variables['mask_rho'][:]
msk_romsv = grd.variables['mask_v'][:]
h_roms = grd.variables['h'][:]
sh2 = msk_roms.shape
etamax, ximax = sh2

scoord = s_coordinate_4(h_roms, theta_b, theta_s, tcline, klevels)
zr = scoord.z_r[:]



a,b=np.meshgrid(x_roms[0,:],zr[:,0,0])

val=np.squeeze(tempini[35,:,0,  :])   #south

fig, ax1 = plt.subplots()
AB=ax1.pcolor(a, zr[:,0,:],val, vmin=val.min(), vmax=val.max(), cmap = 'seismic')
#AB=ax1.pcolor(a, zr[:,0,:], temp[0,::])
ax1.axis('tight')
cbar=fig.colorbar(AB, shrink=0.9, extend='both')
cbar.ax.set_ylabel('v(m/s)', rotation=270)
cbar.ax.get_yaxis().labelpad = 20
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=15)
text.set_font_properties(font)
#plt.savefig('bry_south.png')
plt.ylim(-5000,0)
plt.show()


x_roms = grd.variables['lon_rho'][:]
y_roms = grd.variables['lat_rho'][:]

x_roms=x_roms[:-1,:-1]
y_roms=y_roms[:-1,:-1]



tim=35
lay=10
sp=3


up=file['u'][tim,lay,:,:]
up=up[:-1,:]

vp=file['v'][tim,lay,:,:]
vp=vp[:,:-1]

uu=up

vv=vp
val=np.sqrt((uu**2)+(vv**2))

a=plt.pcolor(x_roms,y_roms, np.squeeze(val),vmin=val.min(),vmax=1)
#a=plt.pcolor(x_roms,y_roms, np.squeeze(val),vmin=val.min(),vmax=val.max())
plt.quiver(x_roms[0:-1:sp, 0:-1:sp],y_roms[0:-1:sp, 0:-1:sp],uu[0:-1:sp, 0:-1:sp], vv[0:-1:sp, 0:-1:sp])
plt.text(-49,-27,figdates[tim].strftime("%m%d - %-I %p"), fontsize=12, fontweight='bold',
        bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 5})
plt.colorbar(a)
plt.show()




val=np.squeeze(tempini[7,:,0,  :])

fig, ax1 = plt.subplots()
AB=ax1.pcolor(a, zr[:,0,:],val, vmin=val.min(), vmax=val.max(), cmap = 'seismic')
#AB=ax1.pcolor(a, zr[:,0,:], temp[0,::])
ax1.axis('tight')
cbar=fig.colorbar(AB, shrink=0.9, extend='both')
cbar.ax.set_ylabel('v(m/s)', rotation=270)
cbar.ax.get_yaxis().labelpad = 20
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=15)
text.set_font_properties(font)
#plt.savefig('bry_south.png')
plt.ylim(-5000,0)

(zind,lonind)=np.where(val==val.max())
ax1.plot(a[zind,lonind], zr[zind,0,lonind], '*', color='green')
(zind,lonind)=np.where(val==val.min())
ax1.plot(a[zind,lonind], zr[zind,0,lonind], '*', color='black')


plt.show()






# ZOOMED FIGURE

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset


fig, ax1 = plt.subplots()
AB=ax1.pcolor(a, zr[:,-1,:], temp[0,::], vmin=-0.8, vmax=0.8, cmap = 'seismic')
#AB=ax1.pcolor(a, zr[:,0,:], temp[0,::])
ax1.axis('tight')
cbar=fig.colorbar(AB, shrink=0.9, extend='both')
cbar.ax.set_ylabel('v(m/s)', rotation=270)
cbar.ax.get_yaxis().labelpad = 20
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=15)
text.set_font_properties(font)
#plt.savefig('bry_south.png')
plt.ylim(-5000,0)

axins = zoomed_inset_axes(ax1, 4, loc=3) # zoom-factor: 2.5, location: upper-left

axins.pcolor(a, zr[:,-1,:], temp[0,::], vmin=-0.8, vmax=0.8, cmap = 'seismic')

x1, x2, y1, y2 = -40, -36, -200, 0 # specify the limits

axins.set_xlim(x1, x2) # apply the x-limits
axins.set_ylim(y1, y2) # apply the y-limits

axins.tick_params(axis="y",direction="in", pad=-22)
axins.tick_params(axis="x",direction="in", pad=-15)

plt.yticks(visible=False)
plt.xticks(visible=False)

mark_inset(ax1, axins, loc1=2, loc2=1, fc="none", ec="0.5")

plt.show()


#loc
'upper right'  : 1,
'upper left'   : 2,
'lower left'   : 3,
'lower right'  : 4,
'right'        : 5,
'center left'  : 6,
'center right' : 7,
'lower center' : 8,
'upper center' : 9,
'center'       : 10


plt.show()





a,b=np.meshgrid(x_roms[0,:],z_fm)

V_int2d = np.load('V_int2d.npy')

plt.pcolor(a,b,V_int2d[:,-1,:], vmin=-0.3, vmax=0.3);plt.ylim(-5000,0);plt.colorbar();plt.show()




fig, ax1 = plt.subplots()
AB=ax1.pcolor(temp[-1,::], vmin=-0.3, vmax=0.3)
ax1.axis('tight')
cbar=fig.colorbar(AB, shrink=0.9, extend='both')
cbar.ax.set_ylabel('v($m^{2}$)', rotation=270)
cbar.ax.get_yaxis().labelpad = 20
text = cbar.ax.yaxis.label
font = matplotlib.font_manager.FontProperties(family='times new roman', style='normal', size=15)
text.set_font_properties(font)
#plt.savefig('bry_south.png')


plt.show()



#######ERA FILE

from netCDF4 import Dataset as dat
from matplotlib import dates
import numpy as np
import matplotlib.pyplot as plt

fileu=dat('u_era.nc')

filev=dat('v_era.nc')


uwind=fileu['Uwind'][:]
vwind=filev['Vwind'][:]


ini=dates.datestr2num('2013-01-01 00:00:00')

datas=ini + fileu['wind_time'][:]/(24*60*60)


ind=np.where(datas==dates.datestr2num('2016-08-15 12:00:00'))

dates.num2date(datas[ind])

sp=2
tim=ind
x_roms=fileu['lon']
y_roms=fileu['lat']


val=np.sqrt((uwind[tim,:,:]**2)+(vwind[tim,:,:]**2))
a=plt.pcolor(x_roms,y_roms, np.squeeze(val),vmin=val.min())
plt.colorbar(a)
plt.quiver(x_roms[0:-1:sp, 0:-1:sp],y_roms[0:-1:sp, 0:-1:sp],np.squeeze(uwind[tim,0:-1:sp, 0:-1:sp]), np.squeeze(vwind[tim,0:-1:sp, 0:-1:sp]))
plt.show()
