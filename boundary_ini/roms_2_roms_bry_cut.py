# -*- coding: utf-8 -*-
#

import numpy as np
import matplotlib.pyplot as plt
from sys import exit
#from gsw import distance
#from gsw import f as fcor
from scipy.io import loadmat
from netCDF4 import Dataset
from vgrid import s_coordinate, s_coordinate_2, s_coordinate_4
from scipy.interpolate import griddata, interp1d
#from oceans.ff_tools import strip_mask
#from seawater import pres, svan
from datetime import datetime as dt
from parameters_bry_in import *
import matplotlib.dates as dates
from funcs_interpolation import horiz_interp_3dvar_cut, vert_interp_3dvar_cut, horiz_interp_3dvar
from funcs_interpolation import vert_interp_3dvar_cut_r2r


##--------------
plt.close('all')


SAVE_INI_BRY_FILES = True


synopsis = 'Bry file Azul project'


## Load ROMS grid.
grd = Dataset(fname_grd_son)
x_roms = grd.variables['lon_rho'][:]
y_roms = grd.variables['lat_rho'][:]
msk_roms = grd.variables['mask_rho'][:]
msk_romsv = grd.variables['mask_v'][:]
h_roms = grd.variables['h'][:]
sh2 = msk_roms.shape
etamax, ximax = sh2
angles=grd['angle'][:]

print('processing grid named ' + fname_grd_son)
print('\n')
print('processing coupling files ' + coup_files[0])

###################
sumtimes=0

for i in range(len(coup_files)): 
  sumtimes = int(Dataset(coup_files[i])['u_eastward'][:].shape[0]) + sumtimes
  
ref_time = str(Dataset(coup_files[0])['ocean_time'].units)  
 
file=Dataset(coup_files[0])
x_fm=file['lon_rho'][:]
x_fm=x_fm[::-1,:]

y_fm=file['lat_rho'][:]
y_fm=y_fm[::-1,:]

if not rotated: 
  minlon = x_fm[0,:] - x_roms.min()
  iml = np.where(np.absolute(minlon)==np.absolute(minlon).min())[0][0]
  maxlon = x_fm[0,:] - x_roms.max()
  imxl = np.where(np.absolute(maxlon)==np.absolute(maxlon).min())[0][0]
 
  minlat = y_fm[:,0] - y_roms.min()
  imla = np.where(np.absolute(minlat)==np.absolute(minlat).min())[0][0]
  maxlat = y_fm[:,0] - y_roms.max()
  imxla = np.where(np.absolute(maxlat)==np.absolute(maxlat).min())[0][0]
  
else:
  minlon = x_fm[0,:] - x_roms[:,0].max()
  iml1 = np.where(np.absolute(minlon)==np.absolute(minlon).min())[0][0]
  minlon = x_fm[-1,:] - x_roms[:,0].max()
  iml2 = np.where(np.absolute(minlon)==np.absolute(minlon).min())[0][0]
  iml=min(iml1,iml2)

  maxlon = x_fm[0,:] - x_roms[:,-1].min()
  imxl1 = np.where(np.absolute(maxlon)==np.absolute(maxlon).min())[0][0]
  maxlon = x_fm[-1,:] - x_roms[:,-1].min()
  imxl2 = np.where(np.absolute(maxlon)==np.absolute(maxlon).min())[0][0]
  imxl=max(imxl1,imxl2) 
 
  maxlat = y_fm[:,0] - y_roms[-1,:].min()
  imxla1 = np.where(np.absolute(maxlat)==np.absolute(maxlat).min())[0][0]
  maxlat = y_fm[:,-1] - y_roms[-1,:].min()
  imxla2 = np.where(np.absolute(maxlat)==np.absolute(maxlat).min())[0][0]
  imxla=min(imxla1,imxla2)

  minlat = y_fm[:,0] - y_roms[0,:].max()
  imla1 = np.where(np.absolute(minlat)==np.absolute(minlat).min())[0][0]
  minlat = y_fm[:,-1] - y_roms[0,:].max()
  imla2 = np.where(np.absolute(minlat)==np.absolute(minlat).min())[0][0]
  imla=max(imla1,imla2)


x_fm = x_fm[imxla - lim:imla + lim,iml - lim:imxl + lim]

y_fm = y_fm[imxla - lim:imla + lim,iml - lim:imxl + lim]

plt.plot(x_fm,y_fm, 'red')   #father grid
plt.plot(x_roms,y_roms, 'blue');plt.show()   #son grid. Check if son is completely inside father

numlevels  = int(file['u_eastward'].shape[1])
  
uclm = np.zeros( (sumtimes,) + (numlevels,) + x_fm.shape)
vclm = np.zeros( (sumtimes,) + (numlevels,) + x_fm.shape)
tempclm = np.zeros( (sumtimes,) + (numlevels,) + x_fm.shape)
saltclm = np.zeros( (sumtimes,) + (numlevels,) + x_fm.shape)
sshclm =  np.zeros( (sumtimes,) + x_fm.shape)
bry_time=np.zeros(sumtimes)


counttime=0
for i in range(len(coup_files)): 
  file=Dataset(coup_files[i])
  u=file['u_eastward'][:]
  u=u[:,:,::-1,:]
  u=u[:,:,imxla - lim:imla + lim,iml - lim:imxl + lim]
  v=file['v_northward'][:]
  v=v[:,:,::-1,:]
  v=v[:,:,imxla - lim:imla + lim,iml - lim:imxl + lim]
  temp=file['temp'][:]
  temp=temp[:,:,::-1,:]
  temp=temp[:,:,imxla - lim:imla + lim,iml - lim:imxl + lim]
  salt=file['salt'][:]
  salt=salt[:,:,::-1,:]
  salt=salt[:,:,imxla - lim:imla + lim,iml - lim:imxl + lim]
  ssh=file['zeta'][:]
  ssh=ssh[:,::-1,:]
  ssh=ssh[:,imxla - lim:imla + lim,iml - lim:imxl + lim]
  numtime  = int(file['u_eastward'].shape[0])
  uclm[counttime:counttime+numtime,::] = u
  vclm[counttime:counttime+numtime,::] = v
  tempclm[counttime:counttime+numtime,::] = temp
  saltclm[counttime:counttime+numtime,::] = salt
  sshclm[counttime:counttime+numtime,::] = ssh
  bry_time[counttime:counttime+numtime] = file['ocean_time'][:]
  counttime=counttime+numtime

uclm[uclm>100]=np.nan
vclm[vclm>100]=np.nan
tempclm[tempclm>100]=np.nan
saltclm[saltclm>100]=np.nan
sshclm[sshclm>100]=np.nan

h_roms_p = file['h'][:]
h_roms_p=h_roms_p[::-1,:]
h_roms_p=h_roms_p[imxla - lim:imla + lim,iml - lim:imxl + lim]

theta_b_p= float(file['theta_b'][:])
theta_s_p=  float(file['theta_s'][:])
tcline_p= float(file['Tcline'][:])
Vstretching_p = float(file['Vstretching'][:])
klevels_p = numlevels

ZETA=np.zeros(x_fm.shape)

if Vstretching_p==4:
  scoord = s_coordinate_4(h_roms_p, theta_b_p, theta_s_p, tcline_p, klevels_p, zeta=ZETA)   #zeta is not used in the computation 
elif Vstretching_p==2:
  scoord = s_coordinate_2(h_roms_p, theta_b_p, theta_s_p, tcline_p, klevels_p, zeta=ZETA)
elif Vstretching_p==1:
  scoord = s_coordinate(h_roms_p, theta_b_p, theta_s_p, tcline_p, klevels_p, zeta=ZETA)

z_parent = scoord.z_r[:]

##interp s layers to son's grid

z_parent_i = np.zeros((z_parent.shape[0],) + x_roms.shape)
for i in range(z_parent.shape[0]):
  z_parent_i[i,:] = griddata((x_fm.ravel(),y_fm.ravel()),z_parent[i,:].ravel(), (x_roms.ravel(),y_roms.ravel())).reshape(x_roms.shape)

z_parent =  z_parent_i 

u_int3d=np.zeros((sumtimes,) + (klevels,) + x_roms.shape)
v_int3d=np.zeros((sumtimes,) + (klevels,) + x_roms.shape)
temp_int3d=np.zeros((sumtimes,) + (klevels,) + x_roms.shape)
salt_int3d=np.zeros((sumtimes,) + (klevels,) + x_roms.shape)
zeta=np.zeros((sumtimes,)  + x_roms.shape)

ubar_int3d=np.zeros((sumtimes,)  + (x_roms.shape[0],) + (x_roms.shape[1]-1,))

vbar_int3d=np.zeros((sumtimes,)  + (x_roms.shape[0]-1,) + (x_roms.shape[1],))


dye=np.zeros((sumtimes,) + (klevels,) + x_roms.shape)

if not rotated:
  rotat= y_fm[:,0] - y_roms[-1,:].min()
  nrot = np.where(np.absolute(rotat)==np.absolute(rotat).min())[0][0]
  rotat= y_fm[:,0] - y_roms[0,:].max()
  srot = np.where(np.absolute(rotat)==np.absolute(rotat).min())[0][0]
  rotat = x_fm[0,:] - x_roms[:,-1].min()
  erot = np.where(np.absolute(rotat)==np.absolute(rotat).min())[0][0]
  rotat = x_fm[0,:] - x_roms[:,0].max()
  wrot = np.where(np.absolute(rotat)==np.absolute(rotat).min())[0][0]
else:
  rotat= y_fm[:,0] - y_roms[-1,:].min()
  nrot1 = np.where(np.absolute(rotat)==np.absolute(rotat).min())[0][0]
  rotat= y_fm[:,-1] - y_roms[-1,:].min()
  nrot2 = np.where(np.absolute(rotat)==np.absolute(rotat).min())[0][0]
  nrot=max(nrot1,nrot2)

  rotat= y_fm[:,0] - y_roms[0,:].max()
  srot1 = np.where(np.absolute(rotat)==np.absolute(rotat).min())[0][0]
  rotat= y_fm[:,-1] - y_roms[0,:].max()
  srot2 = np.where(np.absolute(rotat)==np.absolute(rotat).min())[0][0]
  srot=min(srot1,srot2)

  rotat = x_fm[0,:] - x_roms[:,-1].min()
  erot1 = np.where(np.absolute(rotat)==np.absolute(rotat).min())[0][0]
  rotat = x_fm[-1,:] - x_roms[:,-1].min()
  erot2 = np.where(np.absolute(rotat)==np.absolute(rotat).min())[0][0]
  erot=min(erot1,erot2)

  rotat = x_fm[0,:] - x_roms[:,0].max()
  wrot1 = np.where(np.absolute(rotat)==np.absolute(rotat).min())[0][0]
  rotat = x_fm[-1,:] - x_roms[:,0].max()
  wrot2 = np.where(np.absolute(rotat)==np.absolute(rotat).min())[0][0]
  wrot=max(wrot1, wrot2)
   
   
   
for rt in range(sumtimes):
  #load HYCOM DATA
  print('interpolating time '  + str(rt+1))

  temp_fm = tempclm[rt,:]
  salt_fm = saltclm[rt,:]
  u_fm = uclm[rt,:]
  v_fm = vclm[rt,:]
  ssh_fm = sshclm[rt,:]
     
  
  z_fm = z_parent[:,-1,-1]
  
  if True:
    if cutted:
        ## Horizontally interpolate free-surface height.
        print('notcutted')
        print ("Horizontally interpolate free-surface height.")
        points = (x_fm.ravel(),y_fm.ravel())
        print ('')
        interp_points = (x_roms.ravel(),y_roms.ravel())
        zeta_int2d = griddata(points, ssh_fm.ravel(), interp_points, method='linear').reshape(sh2)
        print ('')
        f = np.isnan(zeta_int2d)
        points = (x_roms[~f],y_roms[~f])
        zeta_int2d = griddata(points, zeta_int2d[~f], interp_points, method='nearest').reshape(sh2)
        np.save('ZETA_int2d.npy',zeta_int2d)
        
        ## Interpolate U, V, TEMP and SALT horizontally to the ROMS grid points.
        U_int2d = np.zeros([u_fm.shape[0], x_roms.shape[0], x_roms.shape[1]])
        a = horiz_interp_3dvar_cut(u_fm, x_fm, y_fm, -z_fm, x_roms[-2:,:], y_roms[-2:,:], h_roms, 'N', nrot)
        U_int2d[:,-2:,:] = a
        a = horiz_interp_3dvar_cut(u_fm, x_fm, y_fm, -z_fm, x_roms[0:2,:], y_roms[0:2,:], h_roms, 'S', srot)
        U_int2d[:,0:2,:] = a
        a = horiz_interp_3dvar_cut(u_fm, x_fm, y_fm, -z_fm, x_roms[:,0:2], y_roms[:,0:2], h_roms, 'W', wrot)
        U_int2d[:,:,0:2] = a
        a = horiz_interp_3dvar_cut(u_fm, x_fm, y_fm, -z_fm, x_roms[:,-2:], y_roms[:,-2:], h_roms, 'E', erot)
        U_int2d[:,:,-2:] = a        
        np.save('U_int2d.npy', U_int2d); del U_int2d

        V_int2d = np.zeros([v_fm.shape[0], x_roms.shape[0], x_roms.shape[1]])
        a = horiz_interp_3dvar_cut(v_fm, x_fm, y_fm, -z_fm, x_roms[-2:,:], y_roms[-2:,:], h_roms, 'N', nrot)
        V_int2d[:,-2:,:] = a
        a = horiz_interp_3dvar_cut(v_fm, x_fm, y_fm, -z_fm, x_roms[0:2,:], y_roms[0:2,:], h_roms, 'S', srot)
        V_int2d[:,0:2,:] = a
        a = horiz_interp_3dvar_cut(v_fm, x_fm, y_fm, -z_fm, x_roms[:,0:2], y_roms[:,0:2], h_roms, 'W', wrot)
        V_int2d[:,:,0:2] = a
        a = horiz_interp_3dvar_cut(v_fm, x_fm, y_fm, -z_fm, x_roms[:,-2:], y_roms[:,-2:], h_roms, 'E', erot)
        V_int2d[:,:,-2:] = a        
        np.save('V_int2d.npy', V_int2d); del V_int2d

          
        TEMP_int2d = np.zeros([temp_fm.shape[0], x_roms.shape[0], x_roms.shape[1]])
        a = horiz_interp_3dvar_cut(temp_fm, x_fm, y_fm, -z_fm, x_roms[-2:,:], y_roms[-2:,:], h_roms, 'N', nrot)
        TEMP_int2d[:,-2:,:] = a
        a = horiz_interp_3dvar_cut(temp_fm, x_fm, y_fm, -z_fm, x_roms[0:2,:], y_roms[0:2,:], h_roms, 'S', srot)
        TEMP_int2d[:,0:2,:] = a
        a = horiz_interp_3dvar_cut(temp_fm, x_fm, y_fm, -z_fm, x_roms[:,0:2], y_roms[:,0:2], h_roms, 'W', wrot)
        TEMP_int2d[:,:,0:2] = a
        a = horiz_interp_3dvar_cut(temp_fm, x_fm, y_fm, -z_fm, x_roms[:,-2:], y_roms[:,-2:], h_roms, 'E', erot)
        TEMP_int2d[:,:,-2:] = a        
        np.save('TEMP_int2d.npy', TEMP_int2d); del TEMP_int2d           
           
        
        SALT_int2d = np.zeros([salt_fm.shape[0], x_roms.shape[0], x_roms.shape[1]])
        a = horiz_interp_3dvar_cut(salt_fm, x_fm, y_fm, -z_fm, x_roms[-2:,:], y_roms[-2:,:], h_roms, 'N', nrot)
        SALT_int2d[:,-2:,:] = a
        a = horiz_interp_3dvar_cut(salt_fm, x_fm, y_fm, -z_fm, x_roms[0:2,:], y_roms[0:2,:], h_roms, 'S', srot)
        SALT_int2d[:,0:2,:] = a
        a = horiz_interp_3dvar_cut(salt_fm, x_fm, y_fm, -z_fm, x_roms[:,0:2], y_roms[:,0:2], h_roms, 'W', wrot)
        SALT_int2d[:,:,0:2] = a
        a = horiz_interp_3dvar_cut(salt_fm, x_fm, y_fm, -z_fm, x_roms[:,-2:], y_roms[:,-2:], h_roms, 'E', erot)
        SALT_int2d[:,:,-2:] = a        
        np.save('SALT_int2d.npy', SALT_int2d); del SALT_int2d        
        
        #exit()
    else:
        print('not cutted')
        print ("Horizontally interpolate free-surface height.")
        points = (x_fm.ravel(),y_fm.ravel())
        print ('')
        interp_points = (x_roms.ravel(),y_roms.ravel())
        zeta_int2d = griddata(points, ssh_fm.ravel(), interp_points, method='linear').reshape(sh2)
        print ('')
        f = np.isnan(zeta_int2d)
        points = (x_roms[~f],y_roms[~f])
        zeta_int2d = griddata(points, zeta_int2d[~f], interp_points, method='nearest').reshape(sh2)
        np.save('ZETA_int2d.npy',zeta_int2d)
        ## Interpolate U, V, TEMP and SALT horizontally to the ROMS grid points.
        U_int2d = horiz_interp_3dvar(u_fm, x_fm, y_fm, -z_fm, x_roms, y_roms, h_roms)
        np.save('U_int2d.npy', U_int2d.filled()); del U_int2d

        V_int2d = horiz_interp_3dvar(v_fm, x_fm, y_fm, -z_fm, x_roms, y_roms, h_roms)
        np.save('V_int2d.npy', V_int2d.filled()); del V_int2d

        TEMP_int2d = horiz_interp_3dvar(temp_fm, x_fm, y_fm, -z_fm, x_roms, y_roms, h_roms)
        np.save('TEMP_int2d.npy', TEMP_int2d.filled()); del TEMP_int2d

        SALT_int2d = horiz_interp_3dvar(salt_fm, x_fm, y_fm, -z_fm, x_roms, y_roms, h_roms)
        np.save('SALT_int2d.npy', SALT_int2d.filled()); del SALT_int2d
        #exit() 
        #exit()
        
  print ('passou oii')

  ## Interpolate U, V, TEMP and SALT vertically to the ROMS grid points.
  ZETA = np.load('ZETA_int2d.npy')

  if Vstretching==4:
        scoord = s_coordinate_4(h_roms, theta_b, theta_s, tcline, klevels, zeta=ZETA)   #zeta is not used in the computation 
  elif Vstretching==2:
        scoord = s_coordinate_2(h_roms, theta_b, theta_s, tcline, klevels, zeta=ZETA)
  elif Vstretching==1:
        scoord = s_coordinate(h_roms, theta_b, theta_s, tcline, klevels, zeta=ZETA)

  zr = scoord.z_r[:]
  

  K,M,L = zr.shape ## Number of full-grid (i.e., non-interior) RHO-points in (sigma,eta,xi) directions.
  Mm = M - 1       ## Number of full-grid (i.e., non-interior) V-points in (eta) direction.
  Lm = L - 1       ## Number of full-grid (i.e., non-interior) U-points in (xi) direction.

  temp = np.nan*np.ones((1,K,M,L))
  salt = np.nan*np.ones((1,K,M,L))
  uRHO = np.nan*np.ones((1,K,M,L))
  vRHO = np.nan*np.ones((1,K,M,L))

  if True:
        # U_int3d, V_int3d, TEMP_int3d, SALT_int3d = vert_interp_3dvar_cut_simultaneous(U_int2d, V_int2d, TEMP_int2d, SALT_int2d, z_parent, zr) 
        
        U_int2d = np.load('U_int2d.npy')
        U_int3d = vert_interp_3dvar_cut(U_int2d, z_parent, zr)
        np.savez('int3d_u.npz', U_int3d=U_int3d)
        # exit()
        del U_int2d, U_int3d

        V_int2d = np.load('V_int2d.npy')
        V_int3d = vert_interp_3dvar_cut(V_int2d, z_parent, zr)
        np.savez('int3d_v.npz', V_int3d=V_int3d)
        # exit()
        del V_int2d, V_int3d

        TEMP_int2d = np.load('TEMP_int2d.npy')
        TEMP_int3d = vert_interp_3dvar_cut_r2r(TEMP_int2d, z_parent, zr)
        np.savez('int3d_t.npz', TEMP_int3d=TEMP_int3d)
        # exit()
        del TEMP_int2d, TEMP_int3d

        SALT_int2d = np.load('SALT_int2d.npy')
        SALT_int3d = vert_interp_3dvar_cut_r2r(SALT_int2d, z_parent, zr)
        np.savez('int3d_s.npz', SALT_int3d=SALT_int3d)
        # exit()
        del SALT_int2d, SALT_int3d

        # np.savez('int3d.npz', U_int3d=U_int3d, V_int3d=V_int3d, TEMP_int3d=TEMP_int3d, SALT_int3d=SALT_int3d)
        #exit()

  if True:
        a = np.load('int3d_u.npz')['U_int3d']
        b = np.load('int3d_v.npz')['V_int3d']
        c = np.load('int3d_t.npz')['TEMP_int3d']
        d = np.load('int3d_s.npz')['SALT_int3d']
        e = np.load('ZETA_int2d.npy')
        
        
  u_int3d[rt,:,:,:]=a
  v_int3d[rt,:,:,:]=b
  temp_int3d[rt,:,:,:]=c
  salt_int3d[rt,:,:,:]=d
  zeta[rt,:,:]=e
  
  angless = np.tile(angles, (a.shape[0],1,1))
  cosa=np.cos(angless)
  sina=np.sin(angless)

  u_int3d[rt,:,:,:]=cosa*a + sina*b
  
  v_int3d[rt,:,:,:]=cosa*b - sina*a   
  
  
## Mask out NaNs.
temp_int3d, salt_int3d, u_int3d, v_int3d, zeta = map(np.ma.masked_invalid, (temp_int3d, salt_int3d, u_int3d, v_int3d, zeta))
print ('ii')

########################



print ('checar')

#CALCULATING BAROTROPIC COMPONENT FOLLOWING

#https://earthscience.stackexchange.com/questions/13167/barotropic-component-definition

ubar_int3d_1=np.zeros((n_time,)  + (x_roms.shape[0],) + (x_roms.shape[1],))

vbar_int3d_1=np.zeros((n_time,)  + (x_roms.shape[0],) + (x_roms.shape[1],))

dzr=abs(np.diff(zr, axis=0, append=0))

for rt in range(int(numdays[0])):
  ubar_int3d_1[rt,:,:] = np.sum(u_int3d[rt,:,:,:]*dzr, axis=0) / abs(zr[0,::])
  vbar_int3d_1[rt,:,:] = np.sum(v_int3d[rt,:,:,:]*dzr, axis=0) / abs(zr[0,::])

## Moving velocity variables to their respective U,V-points.
u_int3d = 0.5*(u_int3d[:,:,:,1:]+u_int3d[:,:,:,:-1])
v_int3d = 0.5*(v_int3d[:,:,1:,:]+v_int3d[:,:,:-1,:])

ubar_int3d = 0.5*(ubar_int3d_1[:,:,1:]+ubar_int3d_1[:,:,:-1])
vbar_int3d = 0.5*(vbar_int3d_1[:,1:,:]+vbar_int3d_1[:,:-1,:])

## Fixing bottom level values of 3D variables (IF NEEDED ONLY).
if False:
        temp_int3d[0,0,:], salt_int3d[0,0,:], u_int3d[0,0,:], v_int3d[0,0,:] = temp_int3d[0,1,:], salt_int3d[0,1,:], u_int3d[0,1,:], v_int3d[0,1,:]

## Masking bad values.
temp_int3d,salt_int3d,zeta,u_int3d,v_int3d = map(np.ma.masked_invalid, (temp_int3d,salt_int3d,zeta,u_int3d,v_int3d))

## Vertically-averaged velocities.  OLD
#for rt in range(sumtimes):
#  ubar_int3d[rt,:,:] = u_int3d[rt,:,:,:].mean(axis=0)
#  vbar_int3d[rt,:,:] = v_int3d[rt,:,:,:].mean(axis=0)
  

##===================
## Plotting to check.
##===================


print ("chegou aqui")


plt.close('all')
if False:
        tempvals = np.arange(0., 30., 2.)
        salvals = np.arange(34., 37.5, 0.5)
        for ieta in xrange(0,etamax,5):
                ietap = ieta + 1
                print ("Plotting line %s of %s"%(ietap,etamax))
                fig, ax = plt.subplots()
                fmsk = msk_roms[ieta,:]==1.
                ax.contourf(np.tile(dstsec[:,fmsk], (klevels,1)), zr[:,ieta,fmsk], temp_int3d[0,:,ieta,fmsk].T, 500)
                cc = ax.contour(np.tile(dstsec[:,fmsk], (klevels,1)), zr[:,ieta,fmsk], temp_int3d[0,:,ieta,fmsk].T, tempvals, colors='k')
                try:
                        ax.clabel(cc, manual=False, fmt='%.0f', color='k')
                except:
                        pass
                cc = ax.contour(np.tile(dstsec[:,fmsk], (klevels,1)), zr[:,ieta,fmsk], salt_int3d[0,:,ieta,fmsk].T, salvals, colors='green')
                try:
                        ax.clabel(cc, manual=False, fmt='%.1f', color='green')
                except:
                        pass
                for k in xrange(temp_int3d.shape[0]):
                        ax.plot(dstsec[k,:], zr[k,ieta,:], color='grey', alpha=.5)

                fmt = 'png'
                if 'sm_' in fname_grd:
                        figname = 'ts_sec_after_interp2roms_smoo/' + str(ietap) + '.' + fmt
                elif 're_' in fname_grd:
                        figname = 'ts_sec_after_interp2roms_real/' + str(ietap) + '.' + fmt

#                plt.savefig(figname, format=fmt, bbox_inches='tight', pad_inches=0.1, dpi=100)
                plt.close('all')

plt.close('all')
if False:
        vvals = np.arange(-1., 1., .025)
        for ieta in xrange(0,etamax-1+5,5):
                ietap = ieta + 1
                print ("Plotting line %s of %s"%(ietap,etamax))
                fig, ax = plt.subplots()
                try:
                        fmsk = msk_romsv[ieta,:]==1.
                except IndexError:
                        ieta = etamax - 2
                        fmsk = msk_romsv[ieta,:]==1.
                vm = np.abs(v_int3d[0,:,ieta,fmsk]).max()
                ax.contourf(np.tile(dstsec[:,fmsk], (klevels,1)), zr[:,ieta,fmsk], v_int3d[0,:,ieta,fmsk].T, 500, vmin=-vm, vmax=vm, cmap=plt.cm.hsv)
                cc = ax.contour(np.tile(dstsec[:,fmsk], (klevels,1)), zr[:,ieta,fmsk], v_int3d[0,:,ieta,fmsk].T, vvals, colors='k')
                try:
                        ax.clabel(cc, manual=False, fmt='%.2f', color='k')
                except:
                        pass
                fmt = 'png'
                if 'sm_' in fname_grd:
                        figname = 'vvel_sec_after_interp2roms_smoo/' + str(ietap) + '.' + fmt
                elif 're_' in fname_grd:
                        figname = 'vvel_sec_after_interp2roms_real/' + str(ietap) + '.' + fmt

#                plt.savefig(figname, format=fmt, bbox_inches='tight', pad_inches=0.1, dpi=100)
                plt.close('all')

##------------------------------------
## Saving *_ini.nc and *_bry.nc files.
##------------------------------------

###------------------------------------------------------------------
### Setting the value of the dye concentration at all vertical levels
### of each cell to the local depth.
###------------------------------------------------------------------



plt.close('all')
# Saving initial conditions file.
if SAVE_INI_BRY_FILES:     
        ##########################################################
        ### - Write boundary conditions netcdf file (*_bry.nc). ###
        ##########################################################
        print ("")
        print ("Writing boundary conditions file.")
        
        ## Renaming variables before saving.
        temp, salt, u, v, ubar, vbar = temp_int3d, salt_int3d, u_int3d, v_int3d, ubar_int3d, vbar_int3d

        if Spherical:
            spherical = 'T'
        else:
            spherical = 'F'

        filenamestr = '_bry.nc'
        filetypestr = 'ROMS boundary conditions file'

        #########################################
        print ("")
        print ("Writing boundary conditions file.")

        ##---- Slicing the 2D fields.
        zeta_south = zeta[:,0,:]
        zeta_east = zeta[:,:,-1]
        zeta_north = zeta[:,-1,:]
        zeta_west = zeta[:,:,0]

        ubar_south = ubar[:,0,:]
        ubar_east = ubar[:,:,-1]
        ubar_north = ubar[:,-1,:]
        ubar_west = ubar[:,:,0]
        

        vbar_south = vbar[:,0,:]
        vbar_east = vbar[:,:,-1]
        vbar_north = vbar[:,-1,:]
        vbar_west = vbar[:,:,0]
        

        ##---- Slicing the 3D fields.
        u_south = u[:,:,0,:]
        u_east = u[:,:,:,-1]
        u_north = u[:,:,-1,:]
        u_west = u[:,:,:,0]
        

        v_south = v[:,:,0,:]
        v_east = v[:,:,:,-1]
        v_north = v[:,:,-1,:]
        v_west = v[:,:,:,0]
        

        temp_south = temp[:,:,0,:]
        temp_east = temp[:,:,:,-1]
        temp_north = temp[:,:,-1,:]
        temp_west = temp[:,:,:,0]
        

        salt_south = salt[:,:,0,:]
        salt_east = salt[:,:,:,-1]
        salt_north = salt[:,:,-1,:]
        salt_west = salt[:,:,:,0]


        dye_south = dye[:,:,0,:]
        dye_east = dye[:,:,:,-1]
        dye_north = dye[:,:,-1,:]
        dye_west = dye[:,:,:,0]
        
        now = dt.now()
        
        outfile = run_name_son + filenamestr
        
        ncfile = Dataset(outfile, mode='w', clobber='true', format='NETCDF3_CLASSIC')
        # creating DIMENSIONS.        

        ncfile.createDimension('xi_rho', size=L)
        ncfile.createDimension('xi_u', size=Lm)
        ncfile.createDimension('xi_v', size=L)
        ncfile.createDimension('eta_rho', size=M)
        ncfile.createDimension('eta_u', size=M)
        ncfile.createDimension('eta_v', size=Mm)
        ncfile.createDimension('s_rho', size=int(klevels))
        ncfile.createDimension('s_w', size=int(klevels+1))
        ncfile.createDimension('zeta_time', size=bry_time.size)
        ncfile.createDimension('temp_time', size=bry_time.size)
        ncfile.createDimension('salt_time', size=bry_time.size)
        ncfile.createDimension('dye_time', size=bry_time.size)
        ncfile.createDimension('v2d_time', size=bry_time.size)
        ncfile.createDimension('v3d_time', size=bry_time.size)
        ncfile.createDimension('one', size=1)

        # creating GLOBAL ATTRIBUTES
        setattr(ncfile, 'type', filetypestr)
        setattr(ncfile, 'title', synopsis)
        setattr(ncfile, 'out_file', outfile)
        setattr(ncfile, 'grd_file', fname_grd)
        setattr(ncfile, 'history', str(now))

        #############################################################################
        ## creating VARIABLES, ATTRIBUTES and ASSIGNING VALUES
        #############################################################################
        ## ---------------------------------------------------------------------------
        ncfile.createVariable('spherical', 'i')
        setattr(ncfile.variables['spherical'], 'long_name', 'Grid type logical switch')
        setattr(ncfile.variables['spherical'], 'flag_values', [0,1])
        setattr(ncfile.variables['spherical'], 'flag_meanings', 'Cartesian spherical')
        ncfile.variables['spherical'][:]  = int(Spherical)

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('Vtransform', 'i', dimensions=('one'))
        setattr(ncfile.variables['Vtransform'], 'long_name', 'Vertical terrain-following transformation equation')
        ncfile.variables['Vtransform'][:] = np.int(Vtransform)

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('Vstretching', 'i', dimensions=('one'))
        setattr(ncfile.variables['Vstretching'], 'long_name', 'Vertical terrain-following stretching function')
        ncfile.variables['Vstretching'][:] = np.int(Vstretching)

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('theta_s', 'd', dimensions=('one'))
        setattr(ncfile.variables['theta_s'], 'long_name', 'S-coordinate surface control parameter')
        ncfile.variables['theta_s'][:] = np.array(scoord.theta_s)

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('theta_b', 'd', dimensions=('one'))
        setattr(ncfile.variables['theta_b'], 'long_name', 'S-coordinate bottom control parameter')
        ncfile.variables['theta_b'][:] = np.array(scoord.theta_b)

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('Tcline', 'd', dimensions=('one'))
        setattr(ncfile.variables['Tcline'], 'long_name', 'S-coordinate surface/bottom layer width')
        setattr(ncfile.variables['Tcline'], 'units', 'm')
        ncfile.variables['Tcline'][:] = np.array(scoord.Tcline)

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('hc', 'd', dimensions=('one'))
        setattr(ncfile.variables['hc'], 'long_name', 'S-coordinate parameter, critical depth')
        setattr(ncfile.variables['hc'], 'units', 'm')
        ncfile.variables['hc'][:]  = np.array(scoord.hc)

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('s_rho', 'd', dimensions=('s_rho'))
        setattr(ncfile.variables['s_rho'], 'long_name', 'S-coordinate at RHO-points')
        setattr(ncfile.variables['s_rho'], 'valid_min', -1.)
        setattr(ncfile.variables['s_rho'], 'valid_max', 0.)
        setattr(ncfile.variables['s_rho'], 'positive', 'up')
        setattr(ncfile.variables['s_rho'], 'standard_name', 'ocean_s_coordinate_g1')
        setattr(ncfile.variables['s_rho'], 'formula_terms', 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc')
        ncfile.variables['s_rho'][:]  = scoord.s_rho

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('s_w', 'd', dimensions=('s_w'))
        setattr(ncfile.variables['s_w'], 'long_name', 'S-coordinate at W-points')
        setattr(ncfile.variables['s_w'], 'valid_min', -1.)
        setattr(ncfile.variables['s_w'], 'valid_max', 0.)
        setattr(ncfile.variables['s_w'], 'positive', 'up')
        setattr(ncfile.variables['s_w'], 'standard_name', 'ocean_s_coordinate_g1')
        setattr(ncfile.variables['s_w'], 'formula_terms', 's: s_rho C: Cs_w eta: zeta depth: h depth_c: hc')
        ncfile.variables['s_w'][:]  = scoord.s_w

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('Cs_r', 'd', dimensions=('s_rho'))
        setattr(ncfile.variables['Cs_r'], 'long_name', 'S-coordinate stretching curve at RHO-points')
        setattr(ncfile.variables['Cs_r'], 'valid_min', -1.)
        setattr(ncfile.variables['Cs_r'], 'valid_max', 0.)
        ncfile.variables['Cs_r'][:]  = scoord.Cs_r

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('Cs_w', 'd', dimensions=('s_w'))
        setattr(ncfile.variables['Cs_w'], 'long_name', 'S-coordinate stretching curve at W-points')
        setattr(ncfile.variables['Cs_w'], 'valid_min', -1.)
        setattr(ncfile.variables['Cs_w'], 'valid_max', 0.)
        ncfile.variables['Cs_w'][:]  = scoord.Cs_w

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('h', 'd', dimensions=('eta_rho', 'xi_rho'))
        setattr(ncfile.variables['h'], 'long_name', 'Final bathymetry at RHO-points')
        setattr(ncfile.variables['h'], 'units', 'm')
        setattr(ncfile.variables['h'], 'coordinates', 'lon_rho lat_rho')
        ncfile.variables['h'][:]  = h_roms

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('lon_rho', 'd', dimensions=('eta_rho', 'xi_rho'))
        setattr(ncfile.variables['lon_rho'], 'long_name', 'Longitude of RHO-points')
        setattr(ncfile.variables['lon_rho'], 'units', 'Degrees_east')
        setattr(ncfile.variables['lon_rho'], 'standard_name', 'longitude')
        ncfile.variables['lon_rho'][:]  = grd.variables['lon_rho'][:]

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('lat_rho', 'd', dimensions=('eta_rho', 'xi_rho'))
        setattr(ncfile.variables['lat_rho'], 'long_name', 'Latitude of RHO-points')
        setattr(ncfile.variables['lat_rho'], 'units', 'Degrees_north')
        setattr(ncfile.variables['lat_rho'], 'standard_name', 'latitude')
        ncfile.variables['lat_rho'][:]  = grd.variables['lat_rho'][:]

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('lon_u', 'd', dimensions=('eta_u', 'xi_u'))
        setattr(ncfile.variables['lon_u'], 'long_name', 'Longitude of U-points')
        setattr(ncfile.variables['lon_u'], 'units', 'Degrees_east')
        setattr(ncfile.variables['lon_u'], 'standard_name', 'longitude')
        ncfile.variables['lon_u'][:]  = grd.variables['lon_u'][:]

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('lat_u', 'd', dimensions=('eta_u', 'xi_u'))
        setattr(ncfile.variables['lat_u'], 'long_name', 'Latitude of U-points')
        setattr(ncfile.variables['lat_u'], 'units', 'Degrees_north')
        setattr(ncfile.variables['lat_u'], 'standard_name', 'latitude')
        ncfile.variables['lat_u'][:]  = grd.variables['lat_u'][:]

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('lon_v', 'd', dimensions=('eta_v', 'xi_v'))
        setattr(ncfile.variables['lon_v'], 'long_name', 'Longitude of V-points')
        setattr(ncfile.variables['lon_v'], 'units', 'Degrees_east')
        setattr(ncfile.variables['lon_v'], 'standard_name', 'longitude')
        ncfile.variables['lon_v'][:]  = grd.variables['lon_v'][:]

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('lat_v', 'd', dimensions=('eta_v', 'xi_v'))
        setattr(ncfile.variables['lat_v'], 'long_name', 'Latitude of V-points')
        setattr(ncfile.variables['lat_v'], 'units', 'Degrees_north')
        setattr(ncfile.variables['lat_v'], 'standard_name', 'latitude')
        ncfile.variables['lat_v'][:]  = grd.variables['lat_v'][:]

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('zeta_time', 'd', dimensions=('zeta_time'))
        setattr(ncfile.variables['zeta_time'], 'long_name', 'Time for free-surface boundary condition (since model initialization)')
        setattr(ncfile.variables['zeta_time'], 'units', ref_time)
        setattr(ncfile.variables['zeta_time'], 'calendar', 'gregorian') # Setting the 'calendar' attribute to 'none' tells ROMS that the forcing fields are perpetual. No re-reading or interpolation is done.
        ncfile.variables['zeta_time'][:] = bry_time

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('temp_time', 'd', dimensions=('temp_time'))
        setattr(ncfile.variables['temp_time'], 'long_name', 'Time for temperature boundary condition (since model initialization)')
        setattr(ncfile.variables['temp_time'], 'units', ref_time)
        setattr(ncfile.variables['temp_time'], 'calendar', 'gregorian') # Setting the 'calendar' attribute to 'none' tells ROMS that the forcing fields are perpetual. No re-reading or interpolation is done.
        ncfile.variables['temp_time'][:] = bry_time

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('salt_time', 'd', dimensions=('salt_time'))
        setattr(ncfile.variables['salt_time'], 'long_name', 'Time for salinity boundary condition (since model initialization)')
        setattr(ncfile.variables['salt_time'], 'units', ref_time)
        setattr(ncfile.variables['salt_time'], 'calendar', 'gregorian') # Setting the 'calendar' attribute to 'none' tells ROMS that the forcing fields are perpetual. No re-reading or interpolation is done.
        ncfile.variables['salt_time'][:] = bry_time

        ncfile.createVariable('dye_time', 'd', dimensions=('dye_time'))
        setattr(ncfile.variables['dye_time'], 'long_name', 'Time for passive tracer boundary condition (since model initialization)')
        setattr(ncfile.variables['dye_time'], 'units', ref_time)
        setattr(ncfile.variables['dye_time'], 'calendar', 'gregorian') # Setting the 'calendar' attribute to 'none' tells ROMS that the forcing fields are perpetual. No re-reading or interpolation is done.
        ncfile.variables['dye_time'][:] = bry_time

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('v2d_time', 'd', dimensions=('v2d_time'))
        setattr(ncfile.variables['v2d_time'], 'long_name', 'Time for 2D momentum boundary condition (since model initialization)')
        setattr(ncfile.variables['v2d_time'], 'units', ref_time)
        setattr(ncfile.variables['v2d_time'], 'calendar', 'gregorian') # Setting the 'calendar' attribute to 'none' tells ROMS that the forcing fields are perpetual. No re-reading or interpolation is done.
        ncfile.variables['v2d_time'][:] = bry_time

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('v3d_time', 'd', dimensions=('v3d_time'))
        setattr(ncfile.variables['v3d_time'], 'long_name', 'Time for 3D momentum boundary condition (since model initialization)')
        setattr(ncfile.variables['v3d_time'], 'units', ref_time)
        setattr(ncfile.variables['v3d_time'], 'calendar', 'gregorian') # Setting the 'calendar' attribute to 'none' tells ROMS that the forcing fields are perpetual. No re-reading or interpolation is done.
        ncfile.variables['v3d_time'][:] = bry_time

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('zeta_south', 'd', dimensions=('zeta_time', 'xi_rho'))
        setattr(ncfile.variables['zeta_south'], 'long_name', 'Free-surface southern boundary condition')
        setattr(ncfile.variables['zeta_south'], 'units', 'm')
        setattr(ncfile.variables['zeta_south'], 'time', 'zeta_time')
        ncfile.variables['zeta_south'][:] = zeta_south

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('zeta_east', 'd', dimensions=('zeta_time', 'eta_rho'))
        setattr(ncfile.variables['zeta_east'], 'long_name', 'Free-surface eastern boundary condition')
        setattr(ncfile.variables['zeta_east'], 'units', 'm')
        setattr(ncfile.variables['zeta_south'], 'time', 'zeta_time')
        ncfile.variables['zeta_east'][:] = zeta_east

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('zeta_west', 'd', dimensions=('zeta_time', 'eta_rho'))
        setattr(ncfile.variables['zeta_west'], 'long_name', 'Free-surface eastern boundary condition')
        setattr(ncfile.variables['zeta_west'], 'units', 'm')
        setattr(ncfile.variables['zeta_west'], 'time', 'zeta_time')
        ncfile.variables['zeta_west'][:] = zeta_west

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('zeta_north', 'd', dimensions=('zeta_time', 'xi_rho'))
        setattr(ncfile.variables['zeta_north'], 'long_name', 'Free-surface northern boundary condition')
        setattr(ncfile.variables['zeta_north'], 'units', 'm')
        setattr(ncfile.variables['zeta_south'], 'time', 'zeta_time')
        ncfile.variables['zeta_north'][:] = zeta_north

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('ubar_south', 'd', dimensions=('v2d_time', 'xi_u'))
        setattr(ncfile.variables['ubar_south'], 'long_name', '2D u-momentum southern boundary condition')
        setattr(ncfile.variables['ubar_south'], 'units', 'meter second-1')
        setattr(ncfile.variables['zeta_south'], 'time', 'v2d_time')
        ncfile.variables['ubar_south'][:] = ubar_south

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('ubar_east', 'd', dimensions=('v2d_time', 'eta_u'))
        setattr(ncfile.variables['ubar_east'], 'long_name', '2D u-momentum eastern boundary condition')
        setattr(ncfile.variables['ubar_east'], 'units', 'meter second-1')
        setattr(ncfile.variables['zeta_south'], 'time', 'v2d_time')
        ncfile.variables['ubar_east'][:] = ubar_east
        
        ## ---------------------------------------------------------------------------
        ncfile.createVariable('ubar_west', 'd', dimensions=('v2d_time', 'eta_u'))
        setattr(ncfile.variables['ubar_west'], 'long_name', '2D u-momentum eastern boundary condition')
        setattr(ncfile.variables['ubar_west'], 'units', 'meter second-1')
        setattr(ncfile.variables['zeta_west'], 'time', 'v2d_time')
        ncfile.variables['ubar_west'][:] = ubar_west

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('ubar_north', 'd', dimensions=('v2d_time', 'xi_u'))
        setattr(ncfile.variables['ubar_north'], 'long_name', '2D u-momentum northern boundary condition')
        setattr(ncfile.variables['ubar_north'], 'units', 'meter second-1')
        setattr(ncfile.variables['zeta_south'], 'time', 'v2d_time')
        ncfile.variables['ubar_north'][:] = ubar_north

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('vbar_south', 'd', dimensions=('v2d_time', 'xi_v'))
        setattr(ncfile.variables['vbar_south'], 'long_name', '2D v-momentum southern boundary condition')
        setattr(ncfile.variables['vbar_south'], 'units', 'meter second-1')
        setattr(ncfile.variables['zeta_south'], 'time', 'v2d_time')
        ncfile.variables['vbar_south'][:] = vbar_south

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('vbar_east', 'd', dimensions=('v2d_time', 'eta_v'))
        setattr(ncfile.variables['vbar_east'], 'long_name', '2D v-momentum eastern boundary condition')
        setattr(ncfile.variables['vbar_east'], 'units', 'meter second-1')
        setattr(ncfile.variables['zeta_south'], 'time', 'v2d_time')
        ncfile.variables['vbar_east'][:] = vbar_east

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('vbar_west', 'd', dimensions=('v2d_time', 'eta_v'))
        setattr(ncfile.variables['vbar_west'], 'long_name', '2D v-momentum eastern boundary condition')
        setattr(ncfile.variables['vbar_west'], 'units', 'meter second-1')
        setattr(ncfile.variables['zeta_west'], 'time', 'v2d_time')
        ncfile.variables['vbar_west'][:] = vbar_west

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('vbar_north', 'd', dimensions=('v2d_time', 'xi_v'))
        setattr(ncfile.variables['vbar_north'], 'long_name', '2D v-momentum northern boundary condition')
        setattr(ncfile.variables['vbar_north'], 'units', 'meter second-1')
        setattr(ncfile.variables['zeta_south'], 'time', 'v2d_time')
        ncfile.variables['vbar_north'][:] = vbar_north

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('temp_south', 'd', dimensions=('temp_time', 's_rho', 'xi_rho'))
        setattr(ncfile.variables['temp_south'], 'long_name', 'Potential temperature southern boundary condition')
        setattr(ncfile.variables['temp_south'], 'units', 'Degrees Celsius')
        setattr(ncfile.variables['temp_south'], 'time', 'temp_time')
        ncfile.variables['temp_south'][:] = temp_south

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('temp_east', 'd', dimensions=('temp_time', 's_rho', 'eta_rho'))
        setattr(ncfile.variables['temp_east'], 'long_name', 'Potential temperature eastern boundary condition')
        setattr(ncfile.variables['temp_east'], 'units', 'Degrees Celsius')
        setattr(ncfile.variables['temp_east'], 'time', 'temp_time')
        ncfile.variables['temp_east'][:] = temp_east

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('temp_west', 'd', dimensions=('temp_time', 's_rho', 'eta_rho'))
        setattr(ncfile.variables['temp_west'], 'long_name', 'Potential temperature eastern boundary condition')
        setattr(ncfile.variables['temp_west'], 'units', 'Degrees Celsius')
        setattr(ncfile.variables['temp_west'], 'time', 'temp_time')
        ncfile.variables['temp_west'][:] = temp_west

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('temp_north', 'd', dimensions=('temp_time', 's_rho', 'xi_rho'))
        setattr(ncfile.variables['temp_north'], 'long_name', 'Potential temperature northern boundary condition')
        setattr(ncfile.variables['temp_north'], 'units', 'Degrees Celsius')
        setattr(ncfile.variables['temp_north'], 'time', 'temp_time')
        ncfile.variables['temp_north'][:] = temp_north

        ncfile.createVariable('dye_south_01', 'd', dimensions=('dye_time', 's_rho', 'xi_rho'))
        setattr(ncfile.variables['dye_south_01'], 'long_name', 'Dye_01 southern boundary condition')
        setattr(ncfile.variables['dye_south_01'], 'units', 'kilogram meter-3')
        setattr(ncfile.variables['dye_south_01'], 'time', 'dye_time')
        ncfile.variables['dye_south_01'][:]  = dye_south

        ncfile.createVariable('dye_east_01', 'd', dimensions=('dye_time', 's_rho', 'eta_rho'))
        setattr(ncfile.variables['dye_east_01'], 'long_name', 'Dye_01 eastern boundary condition')
        setattr(ncfile.variables['dye_east_01'], 'units', 'kilogram meter-3')
        setattr(ncfile.variables['dye_east_01'], 'time', 'dye_time')
        ncfile.variables['dye_east_01'][:]  = dye_east

        ncfile.createVariable('dye_west_01', 'd', dimensions=('dye_time', 's_rho', 'eta_rho'))
        setattr(ncfile.variables['dye_west_01'], 'long_name', 'Dye_01 eastern boundary condition')
        setattr(ncfile.variables['dye_west_01'], 'units', 'kilogram meter-3')
        setattr(ncfile.variables['dye_west_01'], 'time', 'dye_time')
        ncfile.variables['dye_west_01'][:]  = dye_west

        ncfile.createVariable('dye_north_01', 'd', dimensions=('dye_time', 's_rho', 'xi_rho'))
        setattr(ncfile.variables['dye_north_01'], 'long_name', 'Dye_01 northern boundary condition')
        setattr(ncfile.variables['dye_north_01'], 'units', 'kilogram meter-3')
        setattr(ncfile.variables['dye_north_01'], 'time', 'dye_time')
        ncfile.variables['dye_north_01'][:]  = dye_north

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('salt_south', 'd', dimensions=('salt_time', 's_rho', 'xi_rho'))
        setattr(ncfile.variables['salt_south'], 'long_name', 'Practical salinity southern boundary condition')
        setattr(ncfile.variables['salt_south'], 'units', 'Practical salinity units, psu (PSS-78)')
        setattr(ncfile.variables['salt_south'], 'time', 'salt_time')
        ncfile.variables['salt_south'][:] = salt_south

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('salt_east', 'd', dimensions=('salt_time', 's_rho', 'eta_rho'))
        setattr(ncfile.variables['salt_east'], 'long_name', 'Practical salinity eastern boundary condition')
        setattr(ncfile.variables['salt_east'], 'units', 'Practical salinity units, psu (PSS-78)')
        setattr(ncfile.variables['salt_east'], 'time', 'salt_time')
        ncfile.variables['salt_east'][:] = salt_east
        
        ## ---------------------------------------------------------------------------
        ncfile.createVariable('salt_west', 'd', dimensions=('salt_time', 's_rho', 'eta_rho'))
        setattr(ncfile.variables['salt_west'], 'long_name', 'Practical salinity eastern boundary condition')
        setattr(ncfile.variables['salt_west'], 'units', 'Practical salinity units, psu (PSS-78)')
        setattr(ncfile.variables['salt_west'], 'time', 'salt_time')
        ncfile.variables['salt_west'][:] = salt_west
        

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('salt_north', 'd', dimensions=('salt_time', 's_rho', 'xi_rho'))
        setattr(ncfile.variables['salt_north'], 'long_name', 'Practical salinity northern boundary condition')
        setattr(ncfile.variables['salt_north'], 'units', 'Practical salinity units, psu (PSS-78)')
        setattr(ncfile.variables['salt_north'], 'time', 'salt_time')
        ncfile.variables['salt_north'][:] = salt_north

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('u_south', 'd', dimensions=('v3d_time', 's_rho', 'xi_u'))
        setattr(ncfile.variables['u_south'], 'long_name', '3D u-momentum southern boundary condition')
        setattr(ncfile.variables['u_south'], 'units', 'meter second-1')
        setattr(ncfile.variables['u_south'], 'time', 'v3d_time')
        ncfile.variables['u_south'][:] = u_south

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('u_east', 'd', dimensions=('v3d_time', 's_rho', 'eta_u'))
        setattr(ncfile.variables['u_east'], 'long_name', '3D u-momentum eastern boundary condition')
        setattr(ncfile.variables['u_east'], 'units', 'meter second-1')
        setattr(ncfile.variables['u_east'], 'time', 'v3d_time')
        ncfile.variables['u_east'][:] = u_east
        
        ## ---------------------------------------------------------------------------
        ncfile.createVariable('u_west', 'd', dimensions=('v3d_time', 's_rho', 'eta_u'))
        setattr(ncfile.variables['u_west'], 'long_name', '3D u-momentum eastern boundary condition')
        setattr(ncfile.variables['u_west'], 'units', 'meter second-1')
        setattr(ncfile.variables['u_west'], 'time', 'v3d_time')
        ncfile.variables['u_west'][:] = u_west
        

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('u_north', 'd', dimensions=('v3d_time', 's_rho', 'xi_u'))
        setattr(ncfile.variables['u_north'], 'long_name', '3D u-momentum northern boundary condition')
        setattr(ncfile.variables['u_north'], 'units', 'meter second-1')
        setattr(ncfile.variables['u_north'], 'time', 'v3d_time')
        ncfile.variables['u_north'][:] = u_north

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('v_south', 'd', dimensions=('v3d_time', 's_rho', 'xi_v'))
        setattr(ncfile.variables['v_south'], 'long_name', '3D v-momentum southern boundary condition')
        setattr(ncfile.variables['v_south'], 'units', 'meter second-1')
        setattr(ncfile.variables['v_south'], 'time', 'v3d_time')
        ncfile.variables['v_south'][:] = v_south

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('v_east', 'd', dimensions=('v3d_time', 's_rho', 'eta_v'))
        setattr(ncfile.variables['v_east'], 'long_name', '3D v-momentum eastern boundary condition')
        setattr(ncfile.variables['v_east'], 'units', 'meter second-1')
        setattr(ncfile.variables['v_east'], 'time', 'v3d_time')
        ncfile.variables['v_east'][:] = v_east
        
        ## ---------------------------------------------------------------------------
        ncfile.createVariable('v_west', 'd', dimensions=('v3d_time', 's_rho', 'eta_v'))
        setattr(ncfile.variables['v_west'], 'long_name', '3D v-momentum eastern boundary condition')
        setattr(ncfile.variables['v_west'], 'units', 'meter second-1')
        setattr(ncfile.variables['v_west'], 'time', 'v3d_time')
        ncfile.variables['v_west'][:] = v_west
        

        ## ---------------------------------------------------------------------------
        ncfile.createVariable('v_north', 'd', dimensions=('v3d_time', 's_rho', 'xi_v'))
        setattr(ncfile.variables['v_north'], 'long_name', '3D v-momentum northern boundary condition')
        setattr(ncfile.variables['v_north'], 'units', 'meter second-1')
        setattr(ncfile.variables['v_north'], 'time', 'v3d_time')
        ncfile.variables['v_north'][:] = v_north

        #############################################################################
        ncfile.sync()
        ncfile.close()
        print ("Done.")
        print ("")
        
        grd.close()