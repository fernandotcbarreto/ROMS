from scipy.interpolate import griddata, interp1d
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as dat
from vgrid import s_coordinate, s_coordinate_2, s_coordinate_4
import numpy.ma as ma



interp_points = (y_roms.ravel(), x_roms.ravel())  ## Interpolation grid points (the child grid's coordinates).

points = (y_fm.ravel(),  x_fm.ravel())

ndir=griddata(points, v_fm[-1,:,:].ravel(), interp_points, method='linear').reshape(x_roms.shape)

ff=np.isnan(ndir)

interp_points = (y_roms[ff].ravel(), x_roms[ff].ravel())  ## Interpolation grid points (the child grid's coordinates).

points = (y_fm[~ff].ravel(),  x_fm[~ff].ravel())

ndir[ff] = griddata(points, ndir[~ff], interp_points, method='nearest')

#plt.pcolor(x_roms, y_roms,ma.masked_invalid(ndir));plt.colorbar();plt.show()


plt.pcolor(x_roms, y_roms,ndir);plt.colorbar();plt.show()



V_int2d = np.load('V_int2d.npy')
plt.pcolor(x_roms, y_roms,V_int2d[-1,::]);plt.colorbar();plt.show()





interp_points = (x_roms.ravel(), y_roms.ravel())  ## Interpolation grid points (the child grid's coordinates).

points = (x_fm.ravel(),  y_fm.ravel())

ndir=griddata(points, v_fm[-1,:,:].ravel(), interp_points, method='linear').reshape(x_roms.shape)

ff=np.isnan(ndir)

interp_points = (x_roms[ff].ravel(), y_roms[ff].ravel())  ## Interpolation grid points (the child grid's coordinates).

points = (x_fm[~ff].ravel(),  y_fm[~ff].ravel())

ndir[ff] = griddata(points, ndir[~ff], interp_points, method='nearest')

plt.pcolor(x_roms, y_roms,ndir);plt.colorbar();plt.show()




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

a,b=np.meshgrid(x_roms[0,:],zr[:,0,0])
b = np.load('int3d_v.npz')['V_int3d']

plt.pcolor(a, zr[:,-1,:], b[:,-1,:],vmin=-0.3,vmax=0.3)



plt.pcolor(x_roms, y_roms,b[-1,::]);plt.colorbar();plt.show()


plt.pcolor(V_int2d[:,-1,:]);plt.colorbar();plt.show()



plt.pcolor(x_roms, y_roms, ma.masked_invalid(ndir));plt.colorbar();plt.show()





plt.pcolor(x_fm, y_fm,ma.masked_invalid(v_fm[-1,:,:]));plt.colorbar();plt.show()








ma.masked_invalid

@timeit
def loop():
  for i in [-1]:
    ff=np.isnan(v_fm[i,:,:])
    print(ff.mean())
    if ff.mean()==1.0:
      continue
    interp_points = (y_roms.ravel(), x_roms.ravel())  ## Interpolation grid points (the child grid's coordinates).
    points = (y_fm.ravel(),  x_fm.ravel())
    ndir=griddata(points, v_fm[i,:,:].ravel(), interp_points, method='linear').reshape(x_roms.shape)
    ff=np.isnan(ndir)
    interp_points = (y_roms[ff].ravel(), x_roms[ff].ravel())  ## Interpolation grid points (the child grid's coordinates).
    points = (y_roms[~ff].ravel(),  x_roms[~ff].ravel())
    ndir[ff] = griddata(points, ndir[~ff], interp_points, method='nearest')
    print(ndir.shape)
    return ndir
  
a=loop()  

@timeit
def loop1():
  for i in [-1]:
    ff=np.isnan(v_fm[i,:,:])
    print(ff.mean())
    if ff.mean()==1.0:
      continue
    interp_points = (y_roms.ravel(), x_roms.ravel())  ## Interpolation grid points (the child grid's coordinates).
    points = (y_fm[0:10,-10:-1].ravel(),  x_fm[0:10,-10:-1].ravel())
    ndir=griddata(points, v_fm[i,0:10,-10:-1].ravel(), interp_points, method='linear').reshape(x_roms.shape)
    ff=np.isnan(ndir)
    interp_points = (y_roms[ff].ravel(), x_roms[ff].ravel())  ## Interpolation grid points (the child grid's coordinates).
    points = (y_roms[~ff].ravel(),  x_roms[~ff].ravel())
    ndir[ff] = griddata(points, ndir[~ff], interp_points, method='nearest')
    print(ndir.shape)
    return ndir

     
b=loop1()    
    
@timeit
def loop2():
  for i in [-1]:
    ff=np.isnan(v_fm[i,:,:])
    print(ff.mean())
    if ff.mean()==1.0:
      continue
    interp_points = (y_roms.ravel(), x_roms.ravel())  ## Interpolation grid points (the child grid's coordinates).
    points = (y_fm[~ff].ravel(),  x_fm[~ff].ravel())
    ndir=griddata(points, v_fm[-1,~ff].ravel(), interp_points, method='linear').reshape(x_roms.shape)
    ff=np.isnan(ndir)
    interp_points = (y_roms[ff].ravel(), x_roms[ff].ravel())  ## Interpolation grid points (the child grid's coordinates).
    points = (y_roms[~ff].ravel(),  x_roms[~ff].ravel())
    ndir[ff] = griddata(points, ndir[~ff], interp_points, method='nearest')
    print(ndir.shape)
    return ndir
  
c=loop2()  

g=np.where((x_fm[0,:]<-48)|(x_fm[0,:]>-26))

j=np.where((y_fm[:,0]>-17)|(y_fm[:,0]<-26))

x_fm[list(j[0]),list(g[0])].shape

def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        if 'log_time' in kw:
            name = kw.get('log_name', method.__name__.upper())
            kw['log_time'][name] = int((te - ts) * 1000)
        else:
            print '%r  %2.2f ms' % \
                  (method.__name__, (te - ts) * 1000)
        return result
    return timed
    
    
datx=pd.DataFrame(x_roms.copy())
datx.iloc[1:-1,1:-1]=np.nan
x_roms_n=datx.values
x_roms_n=ma.masked_invalid(x_roms_n)

daty=pd.DataFrame(y_roms.copy())
daty.iloc[1:-1,1:-1]=np.nan
y_roms_n=daty.values
y_roms_n = ma.masked_invalid(y_roms_n)



datx=pd.DataFrame(x_fm.copy())     #nao interpola com masked, plano frustado
datx.iloc[8:-8,8:-8]=np.nan
x_fm_n=datx.values
x_fm_n=ma.masked_invalid(x_fm_n)

daty=pd.DataFrame(y_fm.copy())
daty.iloc[8:-8,8:-8]=np.nan
y_fm_n=daty.values
y_fm_n = ma.masked_invalid(y_fm_n)

@timeit
def loop():
  for i in range(v_fm.shape[0]):
    ff=np.isnan(v_fm[i,:,:])
    print(ff.mean())
    if ff.mean()==1.0:
      continue
    interp_points = (y_roms_n.ravel(), x_roms_n.ravel())  ## Interpolation grid points (the child grid's coordinates).
    points = (y_fm_n.ravel(),  x_fm_n.ravel())
    ndir=griddata(points, v_fm[i,:,:].ravel(), interp_points, method='linear').reshape(x_roms.shape)
    ff=np.isnan(ndir)
    interp_points = (y_roms_n[ff].ravel(), x_roms_n[ff].ravel())  ## Interpolation grid points (the child grid's coordinates).
    points = (y_fm_n[~ff].ravel(),  x_fm_n[~ff].ravel())
    ndir[ff] = griddata(points, ndir[~ff], interp_points, method='nearest')