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
from funcs_interpolation import horiz_interp_3dvar, vert_interp_3dvar


##--------------
plt.close('all')


SAVE_INI_BRY_FILES = True

synopsis = 'Ini file Azul project'

## Load ROMS grid.
grd = Dataset(fname_grd)
x_roms = grd.variables['lon_rho'][:]
y_roms = grd.variables['lat_rho'][:]
msk_roms = grd.variables['mask_rho'][:]
msk_romsv = grd.variables['mask_v'][:]
h_roms = grd.variables['h'][:]
sh2 = msk_roms.shape
etamax, ximax = sh2

angle=grd.variables['angle'][:]

cosa=np.cos(angle)
sina=np.sin(angle)

cosa=np.tile(cosa, (klevels,1,1))

sina=np.tile(sina, (klevels,1,1))

############  teste with roms bathymetry

#h_roms[h_roms>6000]=np.nan

#ff=np.isnan(h_roms)

#h_roms=griddata((x_roms[~ff],y_roms[~ff]), h_roms[~ff], (x_roms.ravel(), y_roms.ravel()), method='nearest').reshape(x_roms.shape)

###################

a=[timeini]

numdays= len(a)

n_time=int(numdays)

a = dates.datestr2num(a)

time_vec = np.arange(numdays) + a[0]

saida=np.zeros((n_time,) + (klevels,) + x_roms.shape)
u_int3d=np.zeros((n_time,) + (klevels,) + x_roms.shape)
v_int3d=np.zeros((n_time,) + (klevels,) + x_roms.shape)
temp_int3d=np.zeros((n_time,) + (klevels,) + x_roms.shape)
salt_int3d=np.zeros((n_time,) + (klevels,) + x_roms.shape)
zeta=np.zeros((n_time,)  + x_roms.shape)

ubar_int3d=np.zeros((n_time,)  + (x_roms.shape[0],) + (x_roms.shape[1]-1,))

vbar_int3d=np.zeros((n_time,)  + (x_roms.shape[0]-1,) + (x_roms.shape[1],))


dye=np.zeros((n_time,) + (klevels,) + x_roms.shape)
 
for rt in range(int(numdays)):
  #load HYCOM DATA

  nc_hycom=input_path + dates.num2date( time_vec[rt]).strftime("%Y%m%d")+'.nc'
  print(nc_hycom)


  file=Dataset(nc_hycom)

  x_fm=file['longitude'][:]

  y_fm=file['latitude'][:]
  y_fm=y_fm[::-1] 

  x_fm,y_fm=np.meshgrid(x_fm,y_fm)
  
  minlon = x_fm[0,:] - x_roms.min()
  iml = np.where(np.absolute(minlon)==np.absolute(minlon).min())[0][0]
  maxlon = x_fm[0,:] - x_roms.max()
  imxl = np.where(np.absolute(maxlon)==np.absolute(maxlon).min())[0][0]
  
  minlat = y_fm[:,0] - y_roms.min()
  imla = np.where(np.absolute(minlat)==np.absolute(minlat).min())[0][0]
  maxlat = y_fm[:,0] - y_roms.max()
  imxla = np.where(np.absolute(maxlat)==np.absolute(maxlat).min())[0][0]
  
  x_fm = x_fm[imxla - lim:imla + lim,iml - lim:imxl + lim]
  
  y_fm = y_fm[imxla - lim:imla + lim,iml - lim:imxl + lim]

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


  temp_fm = temp_fm[:,imxla - lim:imla + lim,iml - lim:imxl + lim]
  salt_fm = salt_fm[:,imxla - lim:imla + lim,iml - lim:imxl + lim]
  u_fm = u_fm[:,imxla - lim:imla + lim,iml - lim:imxl + lim]
  v_fm = v_fm[:,imxla - lim:imla + lim,iml - lim:imxl + lim]
  ssh_fm = ssh_fm[imxla - lim:imla + lim,iml - lim:imxl + lim]
     
  
  
  if True:
        ## Horizontally interpolate free-surface height.
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
        etamax, ximax = h_roms.shape
        z_parent = np.expand_dims(z_fm,1)
        z_parent = np.expand_dims(z_parent,1)
        z_parent = np.tile(z_parent, (1,etamax,ximax))
        # U_int3d, V_int3d, TEMP_int3d, SALT_int3d = vert_interp_3dvar_simultaneous(U_int2d, V_int2d, TEMP_int2d, SALT_int2d, z_parent, zr) 
        
        U_int2d = np.load('U_int2d.npy')
        U_int3d = vert_interp_3dvar(U_int2d, z_parent, zr)
        np.savez('int3d_u.npz', U_int3d=U_int3d)
        # exit()
        del U_int2d, U_int3d

        V_int2d = np.load('V_int2d.npy')
        V_int3d = vert_interp_3dvar(V_int2d, z_parent, zr)
        np.savez('int3d_v.npz', V_int3d=V_int3d)
        # exit()
        del V_int2d, V_int3d

        TEMP_int2d = np.load('TEMP_int2d.npy')
        TEMP_int3d = vert_interp_3dvar(TEMP_int2d, z_parent, zr)
        np.savez('int3d_t.npz', TEMP_int3d=TEMP_int3d)
        # exit()
        del TEMP_int2d, TEMP_int3d

        SALT_int2d = np.load('SALT_int2d.npy')
        SALT_int3d = vert_interp_3dvar(SALT_int2d, z_parent, zr)
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
  
  
  u_int3d[rt,:,:,:]=cosa*a + sina*b
  
  v_int3d[rt,:,:,:]=cosa*b - sina*a
  
    
  
  
## Mask out NaNs.
temp_int3d, salt_int3d, u_int3d, v_int3d, zeta = map(np.ma.masked_invalid, (temp_int3d, salt_int3d, u_int3d, v_int3d, zeta))
print ('ii')

########################



print ('checar')


#https://earthscience.stackexchange.com/questions/13167/barotropic-component-definition

ubar_int3d_1=np.zeros((n_time,)  + (x_roms.shape[0],) + (x_roms.shape[1],))

vbar_int3d_1=np.zeros((n_time,)  + (x_roms.shape[0],) + (x_roms.shape[1],))

dzr=abs(np.diff(zr, axis=0, append=0))

for rt in range(int(numdays)):
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

## Vertically-averaged velocities.
#for rt in range(int(numdays)):
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
        filenamestr = '_ini.nc'
        filetypestr = 'ROMS initial conditions file'

        ## Renaming variables before saving.
        temp, salt, u, v, ubar, vbar = temp_int3d, salt_int3d, u_int3d, v_int3d, ubar_int3d, vbar_int3d
#        u=u*0
#        v=v*0
#        ubar=ubar*0
#        vbar=vbar*0
#        zeta=zeta*0
        
        ##########################################################
        ### - Write initial conditions netcdf file (*_ini.nc). ###
        ##########################################################
        print ("")
        print ("Writing initial conditions file.")

        if Spherical:
            spherical = 'T'
        else:
            spherical = 'F'

        ### NOTE: Setting the initial model time to 0.
        t0 = 0.
        outfile = run_name + filenamestr
        ncfile = Dataset(outfile, mode='w', clobber='true', format='NETCDF3_CLASSIC')

        ncfile.createDimension('xi_rho', size=L)
        ncfile.createDimension('xi_u', size=Lm)
        ncfile.createDimension('xi_v', size=L)
        ncfile.createDimension('eta_rho', size=M)
        ncfile.createDimension('eta_u', size=M)
        ncfile.createDimension('eta_v', size=Mm)
        ncfile.createDimension('s_rho', size=int(klevels))
        ncfile.createDimension('s_w', size=int(klevels+1))
        # ncfile.createDimension('tracer', size=2)
        ncfile.createDimension('time', size=1)
        ncfile.createDimension('one', size=1)

        # creating GLOBAL ATTRIBUTES
        now = dt.now()
        setattr(ncfile, 'type', filetypestr)
        setattr(ncfile, 'title', synopsis)
        setattr(ncfile, 'out_file', outfile)
        setattr(ncfile, 'grd_file', fname_grd)
        setattr(ncfile, 'history', str(now))

        #############################################################################
        # creating VARIABLES, ATTRIBUTES and ASSIGNING VALUES

        # ---------------------------------------------------------------------------
        ncfile.createVariable('spherical', 'i')
        setattr(ncfile.variables['spherical'], 'long_name', 'Grid type logical switch')
        setattr(ncfile.variables['spherical'], 'flag_values', [0,1])
        setattr(ncfile.variables['spherical'], 'flag_meanings', 'Cartesian spherical')
        ncfile.variables['spherical'][:]  = int(Spherical)

        # ---------------------------------------------------------------------------
        ncfile.createVariable('Vtransform', 'i', dimensions=('one'))
        setattr(ncfile.variables['Vtransform'], 'long_name', 'Vertical terrain-following transformation equation')
        ncfile.variables['Vtransform'][:]  = np.int(Vtransform)

        # ---------------------------------------------------------------------------
        ncfile.createVariable('Vstretching', 'i', dimensions=('one'))
        setattr(ncfile.variables['Vstretching'], 'long_name', 'Vertical terrain-following stretching function')
        ncfile.variables['Vstretching'][:]  = np.int(Vstretching)

        # ---------------------------------------------------------------------------
        ncfile.createVariable('theta_s', 'd', dimensions=('one'))
        setattr(ncfile.variables['theta_s'], 'long_name', 'S-coordinate surface control parameter')
        ncfile.variables['theta_s'][:]  = np.array(scoord.theta_s)

        # ---------------------------------------------------------------------------
        ncfile.createVariable('theta_b', 'd', dimensions=('one'))
        setattr(ncfile.variables['theta_b'], 'long_name', 'S-coordinate bottom control parameter')
        ncfile.variables['theta_b'][:]  = np.array(scoord.theta_b)

        # ---------------------------------------------------------------------------
        ncfile.createVariable('Tcline', 'd', dimensions=('one'))
        setattr(ncfile.variables['Tcline'], 'long_name', 'S-coordinate surface/bottom layer width')
        setattr(ncfile.variables['Tcline'], 'units', 'm')
        ncfile.variables['Tcline'][:]  = np.array(scoord.Tcline)

        ###### Let ROMS compute hc internally to avoid NETCDF_CHECK_VAR error.
        # ---------------------------------------------------------------------------
        ncfile.createVariable('hc', 'd', dimensions=('one'))
        setattr(ncfile.variables['hc'], 'long_name', 'S-coordinate parameter, critical depth')
        setattr(ncfile.variables['hc'], 'units', 'm')
        ncfile.variables['hc'][:]  = np.array(scoord.hc)

        # ---------------------------------------------------------------------------
        ncfile.createVariable('s_rho', 'd', dimensions=('s_rho'))
        setattr(ncfile.variables['s_rho'], 'long_name', 'S-coordinate at RHO-points')
        setattr(ncfile.variables['s_rho'], 'valid_min', -1.)
        setattr(ncfile.variables['s_rho'], 'valid_max', 0.)
        setattr(ncfile.variables['s_rho'], 'positive', 'up')
        if np.int(Vtransform)==1:
            setattr(ncfile.variables['s_rho'], 'standard_name', 'ocean_s_coordinate_g1')
        elif np.int(Vtransform)==2:
            setattr(ncfile.variables['s_rho'], 'standard_name', 'ocean_s_coordinate_g2')
        setattr(ncfile.variables['s_rho'], 'formula_terms', 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc')
        ncfile.variables['s_rho'][:]  = scoord.s_rho

        # ---------------------------------------------------------------------------
        ncfile.createVariable('s_w', 'd', dimensions=('s_w'))
        setattr(ncfile.variables['s_w'], 'long_name', 'S-coordinate at W-points')
        setattr(ncfile.variables['s_w'], 'valid_min', -1.)
        setattr(ncfile.variables['s_w'], 'valid_max', 0.)
        setattr(ncfile.variables['s_w'], 'positive', 'up')
        if np.int(Vtransform)==1:
            setattr(ncfile.variables['s_rho'], 'standard_name', 'ocean_s_coordinate_g1')
        elif np.int(Vtransform)==2:
            setattr(ncfile.variables['s_rho'], 'standard_name', 'ocean_s_coordinate_g2')
        setattr(ncfile.variables['s_w'], 'formula_terms', 's: s_rho C: Cs_w eta: zeta depth: h depth_c: hc')
        ncfile.variables['s_w'][:]  = scoord.s_w

        # ---------------------------------------------------------------------------
        ncfile.createVariable('Cs_r', 'd', dimensions=('s_rho'))
        setattr(ncfile.variables['Cs_r'], 'long_name', 'S-coordinate stretching curve at RHO-points')
        setattr(ncfile.variables['Cs_r'], 'valid_min', -1.)
        setattr(ncfile.variables['Cs_r'], 'valid_max', 0.)
        ncfile.variables['Cs_r'][:]  = scoord.Cs_r

        # ---------------------------------------------------------------------------
        ncfile.createVariable('Cs_w', 'd', dimensions=('s_w'))
        setattr(ncfile.variables['Cs_w'], 'long_name', 'S-coordinate stretching curve at W-points')
        setattr(ncfile.variables['Cs_w'], 'valid_min', -1.)
        setattr(ncfile.variables['Cs_w'], 'valid_max', 0.)
        ncfile.variables['Cs_w'][:]  = scoord.Cs_w

        # ---------------------------------------------------------------------------
        ncfile.createVariable('h', 'd', dimensions=('eta_rho', 'xi_rho'))
        setattr(ncfile.variables['h'], 'long_name', 'Final bathymetry at RHO-points')
        setattr(ncfile.variables['h'], 'units', 'm')
        setattr(ncfile.variables['h'], 'coordinates', 'lon_rho lat_rho')
        ncfile.variables['h'][:]  = h_roms

        # ---------------------------------------------------------------------------
        ncfile.createVariable('lon_rho', 'd', dimensions=('eta_rho', 'xi_rho'))
        setattr(ncfile.variables['lon_rho'], 'long_name', 'Longitude of RHO-points')
        setattr(ncfile.variables['lon_rho'], 'units', 'Degrees_east')
        setattr(ncfile.variables['lon_rho'], 'standard_name', 'longitude')
        ncfile.variables['lon_rho'][:]  = grd.variables['lon_rho'][:]

        # ---------------------------------------------------------------------------
        ncfile.createVariable('lat_rho', 'd', dimensions=('eta_rho', 'xi_rho'))
        setattr(ncfile.variables['lat_rho'], 'long_name', 'Latitude of RHO-points')
        setattr(ncfile.variables['lat_rho'], 'units', 'Degrees_north')
        setattr(ncfile.variables['lat_rho'], 'standard_name', 'latitude')
        ncfile.variables['lat_rho'][:]  = grd.variables['lat_rho'][:]

        # ---------------------------------------------------------------------------
        ncfile.createVariable('lon_u', 'd', dimensions=('eta_u', 'xi_u'))
        setattr(ncfile.variables['lon_u'], 'long_name', 'Longitude of U-points')
        setattr(ncfile.variables['lon_u'], 'units', 'Degrees_east')
        setattr(ncfile.variables['lon_u'], 'standard_name', 'longitude')
        ncfile.variables['lon_u'][:]  = grd.variables['lon_u'][:]

        # ---------------------------------------------------------------------------
        ncfile.createVariable('lat_u', 'd', dimensions=('eta_u', 'xi_u'))
        setattr(ncfile.variables['lat_u'], 'long_name', 'Latitude of U-points')
        setattr(ncfile.variables['lat_u'], 'units', 'Degrees_north')
        setattr(ncfile.variables['lat_u'], 'standard_name', 'latitude')
        ncfile.variables['lat_u'][:]  = grd.variables['lat_u'][:]

        # ---------------------------------------------------------------------------
        ncfile.createVariable('lon_v', 'd', dimensions=('eta_v', 'xi_v'))
        setattr(ncfile.variables['lon_v'], 'long_name', 'Longitude of V-points')
        setattr(ncfile.variables['lon_v'], 'units', 'Degrees_east')
        setattr(ncfile.variables['lon_v'], 'standard_name', 'longitude')
        ncfile.variables['lon_v'][:]  = grd.variables['lon_v'][:]

        # ---------------------------------------------------------------------------
        ncfile.createVariable('lat_v', 'd', dimensions=('eta_v', 'xi_v'))
        setattr(ncfile.variables['lat_v'], 'long_name', 'Latitude of V-points')
        setattr(ncfile.variables['lat_v'], 'units', 'Degrees_north')
        setattr(ncfile.variables['lat_v'], 'standard_name', 'latitude')
        ncfile.variables['lat_v'][:]  = grd.variables['lat_v'][:]

        # ---------------------------------------------------------------------------
        ncfile.createVariable('ocean_time', 'd', dimensions=('time'))
        setattr(ncfile.variables['ocean_time'], 'long_name', 'Time since model initialization')
        setattr(ncfile.variables['ocean_time'], 'units', 's')
        # setattr(ncfile.variables['ocean_time'], 'calendar', '360.0 days / year')
        ncfile.variables['ocean_time'][:]  = np.array([t0])

        # ---------------------------------------------------------------------------
        ncfile.createVariable('zeta', 'd', dimensions=('time', 'eta_rho', 'xi_rho'))
        setattr(ncfile.variables['zeta'], 'long_name', 'Free-surface elevation anomaly')
        setattr(ncfile.variables['zeta'], 'units', 'm')
        setattr(ncfile.variables['zeta'], 'coordinates', 'lon_rho lat_rho ocean_time')
        ncfile.variables['zeta'][:]  = zeta[0,:,:]

        # ---------------------------------------------------------------------------
        ncfile.createVariable('ubar', 'd', dimensions=('time', 'eta_u', 'xi_u'))
        setattr(ncfile.variables['ubar'], 'long_name', 'Zonal component of the vertically-integrated total velocity')
        setattr(ncfile.variables['ubar'], 'units', 'm s^{-1}')
        setattr(ncfile.variables['ubar'], 'coordinates', 'lon_u lat_u ocean_time')
        ncfile.variables['ubar'][:]  = ubar[0,:,:]

        # ---------------------------------------------------------------------------
        ncfile.createVariable('vbar', 'd', dimensions=('time', 'eta_v', 'xi_v'))
        setattr(ncfile.variables['vbar'], 'long_name', 'Meridional component of the vertically-integrated total velocity')
        setattr(ncfile.variables['vbar'], 'units', 'm s^{-1}')
        setattr(ncfile.variables['vbar'], 'coordinates', 'lon_v lat_v ocean_time')
        ncfile.variables['vbar'][:]  = vbar[0,:,:]

        # ---------------------------------------------------------------------------
        ncfile.createVariable('u', 'd', dimensions=('time', 's_rho', 'eta_u', 'xi_u'))
        setattr(ncfile.variables['u'], 'long_name', 'Zonal component of the total velocity')
        setattr(ncfile.variables['u'], 'units', 'm s^{-1}')
        setattr(ncfile.variables['u'], 'coordinates', 'lon_u lat_u s_rho ocean_time')
        ncfile.variables['u'][:]  = u[0,:,:,:]

        # ---------------------------------------------------------------------------
        ncfile.createVariable('v', 'd', dimensions=('time', 's_rho', 'eta_v', 'xi_v'))
        setattr(ncfile.variables['v'], 'long_name', 'Meridional component of the total velocity')
        setattr(ncfile.variables['v'], 'units', 'm s^{-1}')
        setattr(ncfile.variables['v'], 'coordinates', 'lon_v lat_v s_rho ocean_time')
        ncfile.variables['v'][:]  = v[0,:,:,:]

        # ---------------------------------------------------------------------------
        ncfile.createVariable('temp', 'd', dimensions=('time', 's_rho', 'eta_rho', 'xi_rho'))
        setattr(ncfile.variables['temp'], 'long_name', 'Potential temperature')
        setattr(ncfile.variables['temp'], 'units', 'Degrees Celsius')
        setattr(ncfile.variables['temp'], 'coordinates', 'lon_rho lat_rho s_rho ocean_time')
        ncfile.variables['temp'][:]  = temp[0,:,:,:]

        # ---------------------------------------------------------------------------
        ncfile.createVariable('salt', 'd', dimensions=('time', 's_rho', 'eta_rho', 'xi_rho'))
        setattr(ncfile.variables['salt'], 'long_name', 'Practical salinity')
        setattr(ncfile.variables['salt'], 'units', 'Practical salinity units, psu (PSS-78)')
        setattr(ncfile.variables['salt'], 'coordinates', 'lon_rho lat_rho s_rho ocean_time')
        ncfile.variables['salt'][:]  = salt[0,:,:,:]

        # ---------------------------------------------------------------------------
        ncfile.createVariable('dye_01', 'd', dimensions=('time', 's_rho', 'eta_rho', 'xi_rho'))
        setattr(ncfile.variables['dye_01'], 'long_name', 'Dye concentration')
        setattr(ncfile.variables['dye_01'], 'units', 'kilogram meter-3')
        setattr(ncfile.variables['dye_01'], 'coordinates', 'lon_rho lat_rho s_rho ocean_time')
        ncfile.variables['dye_01'][:]  = dye[0,:,:,:]

        #############################################################################
        ncfile.sync()
        ncfile.close()
        print ("Done.")
        print ("")

        grd.close()