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


def extrap1d(interpolator):
    """
    http://stackoverflow.com/questions/2745329/
    How to make scipy.interpolate give an extrapolated result beyond the
    input range.
    """
    xs, ys = interpolator.x, interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0] + (x - xs[0]) * (ys[1] - ys[0]) / (xs[1] - xs[0])
        elif x > xs[-1]:
            return (ys[-1] + (x - xs[-1]) * (ys[-1] - ys[-2]) /
                    (xs[-1] - xs[-2]))
        else:
            return interpolator(x)

    def ufunclike(xs):
        return np.array(map(pointwise, np.array(xs)))

    return ufunclike


def horiz_interp_3dvar(variable_pgrd, lon_pgrd, lat_pgrd, z_pgrd, lon_cgrd, lat_cgrd, topo_cgrd):
	print ("Horizontally interpolating 3D variable...")

	interp_points = (lat_cgrd.ravel(), lon_cgrd.ravel())  ## Interpolation grid points (the child grid's coordinates).

	N = variable_pgrd.shape[0]
	Nn = np.sum(z_pgrd<=topo_cgrd.max())
	####
	
	
	Nn=N
	###
	points = (lat_pgrd.ravel(), lon_pgrd.ravel())	

	variable_cgrd = np.ma.zeros( (Nn,)+topo_cgrd.shape )
	
	for k in reversed(range(Nn)):  
           pvar3d_slice = variable_pgrd[k,:,:]
           ff=np.isnan(pvar3d_slice)
           if ff.mean()==1.0:
              ndir=variable_cgrd[k+1,:,:]
           else:    
              try:
                points = (lat_pgrd[~ff], lon_pgrd[~ff])
                interp_points = (lat_cgrd.ravel(), lon_cgrd.ravel())
                ndir=griddata(points, pvar3d_slice[~ff], interp_points, method='linear').reshape(lon_cgrd.shape)
              except:
                points = (lat_pgrd[~ff], lon_pgrd[~ff])
                interp_points = (lat_cgrd.ravel(), lon_cgrd.ravel())
                ndir=griddata(points, pvar3d_slice[~ff], interp_points, method='nearest').reshape(lon_cgrd.shape)
                
              ff=np.isnan(ndir)
              if ff.mean()==1.0:
                 ndir=variable_cgrd[k+1,:,:]   
                 
              ff=np.isnan(ndir)   
              if np.sum(ff)!= 0:
                  points=(lat_cgrd[~ff], lon_cgrd[~ff])	
                  interp_points = (lat_cgrd[ff], lon_cgrd[ff])   
                  ndir[ff] = griddata(points, ndir[~ff], interp_points, method='nearest')
	          
           
           variable_cgrd[k,:,:] = ndir  
                   
	return variable_cgrd


def horiz_interp_3dvar_cut(variable_pgrd, lon_pgrd, lat_pgrd, z_pgrd, lon_cgrd, lat_cgrd, topo_cgrd, sec, rotlim):
	print ("Horizontally interpolating 3D variable...")
        

	N = variable_pgrd.shape[0]
	####
	
	
	Nn=N
	###

	variable_cgrd = np.ma.zeros( (Nn,)+lon_cgrd.shape )
		
	for k in reversed(range(Nn)):  
              pvar3d_slice = variable_pgrd[k,:,:]
           
              if sec is 'N':
                blk = np.isnan(pvar3d_slice[:lim*2 + rotlim,:])
                if blk.mean()==1.0:
                  if (k==Nn-1):
                    blk = np.isnan(pvar3d_slice[:,:])              
                    points = (lat_pgrd[~blk], lon_pgrd[~blk])
                    interp_points = (lat_cgrd.ravel(), lon_cgrd.ravel())
                    ndir=griddata(points, pvar3d_slice[~blk], interp_points, method='nearest').reshape(lon_cgrd.shape)                 
                  else:
                    ndir=variable_cgrd[k+1,:,:]
                else:     
                  try:
                    points = (lat_pgrd[:lim*2 + rotlim,:][~blk], lon_pgrd[:lim*2+ rotlim,:][~blk])
                    interp_points = (lat_cgrd.ravel(), lon_cgrd.ravel())
                    ndir=griddata(points, pvar3d_slice[:lim*2+ rotlim,:][~blk], interp_points, method='linear').reshape(lon_cgrd.shape)
                  except:
                    points = (lat_pgrd[:lim*2+ rotlim,:][~blk], lon_pgrd[:lim*2+ rotlim,:][~blk])
                    interp_points = (lat_cgrd.ravel(), lon_cgrd.ravel())
                    ndir=griddata(points, pvar3d_slice[:lim*2+ rotlim,:][~blk], interp_points, method='nearest').reshape(lon_cgrd.shape)        
                    
                ff=np.isnan(ndir)
                if ff.mean()==1.0:
                   try:
                     ndir=variable_cgrd[k+1,:,:]   
                   except:
                     blk = np.isnan(pvar3d_slice[:,:])              
                     points = (lat_pgrd[~blk], lon_pgrd[~blk])
                     interp_points = (lat_cgrd.ravel(), lon_cgrd.ravel())
                     ndir=griddata(points, pvar3d_slice[~blk], interp_points, method='nearest').reshape(lon_cgrd.shape)
                    
                ff=np.isnan(ndir)   
                if np.sum(ff)!= 0:
                   points=(lat_cgrd[~ff], lon_cgrd[~ff])	
                   interp_points = (lat_cgrd[ff], lon_cgrd[ff])   
                   ndir[ff] = griddata(points, ndir[~ff], interp_points, method='nearest')
	           
              if sec is 'S':
                blk = np.isnan(pvar3d_slice[-lim*2+ rotlim:,:])              
                if blk.mean()==1.0:
                  if (k==Nn-1):
                    blk = np.isnan(pvar3d_slice[:,:])              
                    points = (lat_pgrd[~blk], lon_pgrd[~blk])
                    interp_points = (lat_cgrd.ravel(), lon_cgrd.ravel())
                    ndir=griddata(points, pvar3d_slice[~blk], interp_points, method='nearest').reshape(lon_cgrd.shape)                 
                  else:
                    ndir=variable_cgrd[k+1,:,:]
                else:
                  try:
                    points = (lat_pgrd[-lim*2+ rotlim:,:][~blk], lon_pgrd[-lim*2+ rotlim:,:][~blk])
                    interp_points = (lat_cgrd.ravel(), lon_cgrd.ravel())
                    ndir=griddata(points, pvar3d_slice[-lim*2+ rotlim:,:][~blk], interp_points, method='linear').reshape(lon_cgrd.shape)
                  except:
                    points = (lat_pgrd[-lim*2+ rotlim:,:][~blk], lon_pgrd[-lim*2+ rotlim:,:][~blk])
                    interp_points = (lat_cgrd.ravel(), lon_cgrd.ravel())
                    ndir=griddata(points, pvar3d_slice[-lim*2+ rotlim:,:][~blk], interp_points, method='nearest').reshape(lon_cgrd.shape)
              
                ff=np.isnan(ndir)
                if ff.mean()==1.0:
                   try:
                     ndir=variable_cgrd[k+1,:,:]   
                   except:
                     blk = np.isnan(pvar3d_slice[:,:])              
                     points = (lat_pgrd[~blk], lon_pgrd[~blk])
                     interp_points = (lat_cgrd.ravel(), lon_cgrd.ravel())
                     ndir=griddata(points, pvar3d_slice[~blk], interp_points, method='nearest').reshape(lon_cgrd.shape)
                     
                ff=np.isnan(ndir)   
                if np.sum(ff)!= 0:
                   points=(lat_cgrd[~ff], lon_cgrd[~ff])	
                   interp_points = (lat_cgrd[ff], lon_cgrd[ff])   
                   ndir[ff] = griddata(points, ndir[~ff], interp_points, method='nearest')
	           
              if sec is 'W':
                blk = np.isnan(pvar3d_slice[:,:lim*2+ rotlim])  
                if blk.mean()==1.0:
                  if (k==Nn-1):
                    blk = np.isnan(pvar3d_slice[:,:])              
                    points = (lat_pgrd[~blk], lon_pgrd[~blk])
                    interp_points = (lat_cgrd.ravel(), lon_cgrd.ravel())
                    ndir=griddata(points, pvar3d_slice[~blk], interp_points, method='nearest').reshape(lon_cgrd.shape)                 
                  else:
                    ndir=variable_cgrd[k+1,:,:]
                else:
                  try:
                    points = (lat_pgrd[:,:lim*2+ rotlim][~blk], lon_pgrd[:,:lim*2+ rotlim][~blk])
                    interp_points = (lat_cgrd.ravel(), lon_cgrd.ravel())
                    ndir=griddata(points, pvar3d_slice[:,:lim*2+ rotlim][~blk], interp_points, method='linear').reshape(lon_cgrd.shape)
                  except:
                    points = (lat_pgrd[:,:lim*2+ rotlim][~blk], lon_pgrd[:,:lim*2+ rotlim][~blk])
                    interp_points = (lat_cgrd.ravel(), lon_cgrd.ravel())                 
                    ndir=griddata(points, pvar3d_slice[:,:lim*2+ rotlim][~blk], interp_points, method='nearest').reshape(lon_cgrd.shape)                    
              
                ff=np.isnan(ndir)
                if ff.mean()==1.0:
                   try:
                     ndir=variable_cgrd[k+1,:,:]   
                   except:
                     blk = np.isnan(pvar3d_slice[:,:])              
                     points = (lat_pgrd[~blk], lon_pgrd[~blk])
                     interp_points = (lat_cgrd.ravel(), lon_cgrd.ravel())
                     ndir=griddata(points, pvar3d_slice[~blk], interp_points, method='nearest').reshape(lon_cgrd.shape)                     
                 
                ff=np.isnan(ndir)   
                if np.sum(ff)!= 0:
                   points=(lat_cgrd[~ff], lon_cgrd[~ff])	
                   interp_points = (lat_cgrd[ff], lon_cgrd[ff])   
                   ndir[ff] = griddata(points, ndir[~ff], interp_points, method='nearest')
	           
              if sec is 'E':
                blk = np.isnan(pvar3d_slice[:,-lim*2+ rotlim:])  
                if blk.mean()==1.0:
                  if (k==Nn-1):
                    blk = np.isnan(pvar3d_slice[:,:])              
                    points = (lat_pgrd[~blk], lon_pgrd[~blk])
                    interp_points = (lat_cgrd.ravel(), lon_cgrd.ravel())
                    ndir=griddata(points, pvar3d_slice[~blk], interp_points, method='nearest').reshape(lon_cgrd.shape)                 
                  else:
                    ndir=variable_cgrd[k+1,:,:]
                else: 
                  try:
                    points = (lat_pgrd[:,-lim*2+ rotlim:][~blk], lon_pgrd[:,-lim*2+ rotlim:][~blk])
                    interp_points = (lat_cgrd.ravel(), lon_cgrd.ravel())
                    ndir=griddata(points, pvar3d_slice[:,-lim*2+ rotlim:][~blk], interp_points, method='linear').reshape(lon_cgrd.shape)
                  except:
                    points = (lat_pgrd[:,-lim*2+ rotlim:][~blk], lon_pgrd[:,-lim*2+ rotlim:][~blk])
                    interp_points = (lat_cgrd.ravel(), lon_cgrd.ravel())
                    ndir=griddata(points, pvar3d_slice[:,-lim*2+ rotlim:][~blk], interp_points, method='nearest').reshape(lon_cgrd.shape)  
                    
                ff=np.isnan(ndir)
                if ff.mean()==1.0:
                   try:
                     ndir=variable_cgrd[k+1,:,:]   
                   except:
                     blk = np.isnan(pvar3d_slice[:,:])              
                     points = (lat_pgrd[~blk], lon_pgrd[~blk])
                     interp_points = (lat_cgrd.ravel(), lon_cgrd.ravel())
                     ndir=griddata(points, pvar3d_slice[~blk], interp_points, method='nearest').reshape(lon_cgrd.shape)
                     
                ff=np.isnan(ndir)   
                if np.sum(ff)!= 0:
                   points=(lat_cgrd[~ff], lon_cgrd[~ff])	
                   interp_points = (lat_cgrd[ff], lon_cgrd[ff])   
                   ndir[ff] = griddata(points, ndir[~ff], interp_points, method='nearest')	           
           
           
              variable_cgrd[k,::] = ndir  
	return variable_cgrd


def vert_interp_3dvar_cut(pvar3d, z_parent, z_child):
	global zp, fp1, zc
	print ("Vertically interpolating 3D variable...")
	cvar3d = np.zeros((z_child.shape[0],) + pvar3d.shape[1:])

	zz_parent = np.unique(z_parent)

	etamax, ximax = pvar3d.shape[1:]
	
	for ieta in [0, 1,etamax-2,etamax-1]:
		for ixi in range(ximax):

		     zp = z_parent[:,ieta,ixi]
		     zc = z_child[:,ieta,ixi]
		     fp = pvar3d[:,ieta,ixi]

		     if noextra:
                        fi = np.interp(zc,zp,fp)  
		     else:
                        
		        I = interp1d(zp, fp, kind='linear', bounds_error=False, assume_sorted=False)		
		        I = extrap1d(I)
			
		        if ixi==6 and ieta==6:
                          print ('ZP',fp)
                                                  
		        fi = I(zc)
			
		     cvar3d[:,ieta,ixi] = fi
	             
	for ixi in [0,1,ximax-2, ximax-1]:
	   for ieta in range(etamax):

		     zp = z_parent[:,ieta,ixi]
		     zc = z_child[:,ieta,ixi]
		     fp = pvar3d[:,ieta,ixi]
			
		     if noextra:
                        fi = np.interp(zc,zp,fp)
		     else:
                        
		        I = interp1d(zp, fp, kind='linear', bounds_error=False, assume_sorted=False)		
		        I = extrap1d(I)
			
		        if ixi==6 and ieta==6:
                          print ('ZP',fp)
                                                  
		        fi = I(zc)
		     cvar3d[:,ieta,ixi] = fi
        
	return cvar3d





def vert_interp_3dvar(pvar3d, z_parent, z_child):
	global zp, fp1, zc
	print ("Vertically interpolating 3D variable...")
	cvar3d = np.zeros((z_child.shape[0],) + pvar3d.shape[1:])

	zz_parent = np.unique(z_parent)

	etamax, ximax = pvar3d.shape[1:]
	for ieta in range(etamax):
		# print ("")
		np.disp("Row %s of %s"%(str(ieta+1),str(etamax)))
		# print ("")
		for ixi in range(ximax):
			# if ixi%20==0:
			# 	np.disp("Col %s of %s"%(str(ixi+1),str(ximax)))
			# np.disp("Col %s of %s"%(str(ixi+1),str(ximax)))

		     zp = z_parent[:,ieta,ixi]
		     zc = z_child[:,ieta,ixi]
		     fp = pvar3d[:,ieta,ixi]

		     if noextra:
                        fi = np.interp(zc,zp,fp)
		     else:
			# If no parent variable deeper than the first child s-level,
			# Search horizontally on the parent grid for the nearest value.

		        I = interp1d(zp, fp, kind='linear', bounds_error=False, assume_sorted=False)		
		        I = extrap1d(I)
			
		        if ixi==6 and ieta==6:
                          print ('ZP',fp)
                                                  
		        fi = I(zc)
			
		     cvar3d[:,ieta,ixi] = fi

	return cvar3d

def vert_interp_3dvar_r2r(pvar3d, z_parent, z_child):
	global zp, fp1, zc
	print ("Vertically interpolating 3D variable...")
	cvar3d = np.zeros((z_child.shape[0],) + pvar3d.shape[1:])

	zz_parent = np.unique(z_parent)

	etamax, ximax = pvar3d.shape[1:]
	for ieta in range(etamax):
		# print ("")
		np.disp("Row %s of %s"%(str(ieta+1),str(etamax)))
		# print ("")
		for ixi in range(ximax):
			# if ixi%20==0:
			# 	np.disp("Col %s of %s"%(str(ixi+1),str(ximax)))
			# np.disp("Col %s of %s"%(str(ixi+1),str(ximax)))

		     zp = z_parent[:,ieta,ixi]
		     zc = z_child[:,ieta,ixi]
		     fp = pvar3d[:,ieta,ixi]

		     if noextra:
                        fi = np.interp(zc,zp,fp)
#                        if ( (zc.min() < zp.min()) and (abs(zc.min() - zp.min()) > 5) ):
                        if ( (zc.min() < zp.min()) ):
                          for ki in range(len(zc)):
                            if (zc[ki] < z_parent[:,ieta,:].min()):
                              zc[ki] = z_parent[:,ieta,:].min() 
                            if (zc[ki] < zp.min()): 
                              for ixi2 in range(ximax):
                                if (z_parent[:,ieta,ixi2].min() <= zc[ki]):
                                  fi[ki] =  np.interp(zc[ki],z_parent[:,ieta,ixi2],pvar3d[:,ieta,ixi2])
                                  break
                            else:
                              break  
		     else:
			# If no parent variable deeper than the first child s-level,
			# Search horizontally on the parent grid for the nearest value.

		        I = interp1d(zp, fp, kind='linear', bounds_error=False, assume_sorted=False)		
		        I = extrap1d(I)
			
		        if ixi==6 and ieta==6:
                          print ('ZP',fp)
                                                  
		        fi = I(zc)
			
		     cvar3d[:,ieta,ixi] = fi

	return cvar3d


def vert_interp_3dvar_cut_r2r(pvar3d, z_parent, z_child):
	global zp, fp1, zc
	print ("Vertically interpolating 3D variable...")
	cvar3d = np.zeros((z_child.shape[0],) + pvar3d.shape[1:])

	zz_parent = np.unique(z_parent)

	etamax, ximax = pvar3d.shape[1:]
	
	for ieta in [0, 1,etamax-2,etamax-1]:
		for ixi in range(ximax):

		     zp = z_parent[:,ieta,ixi]
		     zc = z_child[:,ieta,ixi]
		     fp = pvar3d[:,ieta,ixi]

		     if noextra:
                        fi = np.interp(zc,zp,fp) 
#                        if ( (zc.min() < zp.min()) and (abs(zc.min() - zp.min()) > 5) ):
                        if ( (zc.min() < zp.min()) ):
                          for ki in range(len(zc)):
                            if (zc[ki] < z_parent[:,ieta,:].min()):
                              zc[ki] = z_parent[:,ieta,:].min() 
                            if (zc[ki] < zp.min()): 
                              for ixi2 in range(ximax):
                                if (z_parent[:,ieta,ixi2].min() <= zc[ki]):
                                  fi[ki] =  np.interp(zc[ki],z_parent[:,ieta,ixi2],pvar3d[:,ieta,ixi2])
                                  break
                            else:
                              break  
		     else:
                        
		        I = interp1d(zp, fp, kind='linear', bounds_error=False, assume_sorted=False)		
		        I = extrap1d(I)
			
		        if ixi==6 and ieta==6:
                          print ('ZP',fp)
                                                  
		        fi = I(zc)
			
		     cvar3d[:,ieta,ixi] = fi
	             
	for ixi in [0,1,ximax-2, ximax-1]:
	   for ieta in range(etamax):

		     zp = z_parent[:,ieta,ixi]
		     zc = z_child[:,ieta,ixi]
		     fp = pvar3d[:,ieta,ixi]
			
		     if noextra:
                        fi = np.interp(zc,zp,fp)
#                        if ( (zc.min() < zp.min()) and (abs(zc.min() - zp.min()) > 5) ):
                        if ( (zc.min() < zp.min()) ):
                          for ki in range(len(zc)):
                            if (zc[ki] < z_parent[:,:,ixi].min()):
                              zc[ki] = z_parent[:,:,ixi].min() 
                            if (zc[ki] < zp.min()): 
                              for ieta2 in range(etamax):
                                if (z_parent[:,ieta2,ixi].min() <= zc[ki]):
                                  fi[ki] =  np.interp(zc[ki],z_parent[:,ieta2,ixi],pvar3d[:,ieta2,ixi])
                                  break
                            else:
                              break  
		     else:
                        
		        I = interp1d(zp, fp, kind='linear', bounds_error=False, assume_sorted=False)		
		        I = extrap1d(I)
			
		        if ixi==6 and ieta==6:
                          print ('ZP',fp)
                                                  
		        fi = I(zc)
		     cvar3d[:,ieta,ixi] = fi
        
	return cvar3d