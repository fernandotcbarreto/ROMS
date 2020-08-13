#
#  D_NUDGCOEF:  Driver script to create a ROMS nudging coefficients file.
#
#  This a user modifiable script that can be used to prepare ROMS
#  nudging inverse time scales NetCDF file.  It sets-up  all the
#  necessary parameters and variables. USERS can use this as a
#  prototype for their application.
#
#  Nudging to climatology can be used in ROMS for various purposes:
#
#  (1) Improve the behavior of open boundary conditions.
#  (2) Used in conjunction with sponges.
#  (3) Minimize numerical diapycnal mixing of tracers over steep
#      bathymetry (improve T-S properties in deep water masses).  For
#      example, we can nudge to T-S climatology is areas depeer than
#      1500 m.
#
#  The inverse nudging coefficients have units of 1/time.  The default
#  input units in ROMS is 1/day but 1/second is also possible. The
#  routine 'get_nudgcoef.F' will check the 'units' attribute to compute
#  the conversion factor for 1/second. Users need to be sure that the
#  'units' variable attribute is consistent with the written data.
#
#  The variable names for the nudging coefficients is as follows:
#
#     M2_NudgeCoef       for 2D momentum
#     M3_NudgeCoef       for 3D momentum
#     temp_NudgeCoef     for potential temperature
#     salt_NudgeCoef     for salinity
#     ...
#     NO3_NudgeCoef      for nitrate
#     ...
#     tracer_NudgeCoef   for any generic tracer
#
#  They are all defined at RHO-points. If the nudging coefficients for
#  a specific tracer are available in the NetCDF, ROMS will read that
#  NetCDF variable. If NOT and the generic coefficients 'tracer_NudgeCoef'
#  are available, ROMS will process those values instead.
#
#  Notice that the input swicth 'LnudgeTCLM(itrc,ng)' in ROMS input
#  script 'ocean.in' will control which tracer to nudge in the desired
#  grid.
#
#  Currently, the nudging coefficients are time invariant in ROMS.  The
#  same scales are used for the entire simulation.
#

# svn $Id: d_nudgcoef.m 796 2016-05-11 01:45:21Z arango $
#=========================================================================#
#  Copyright (c) 2002-2016 The ROMS/TOMS Group                            #
#    Licensed under a MIT/X style license                                 #
#    See License_ROMS.txt                           Hernan G. Arango      #
#=========================================================================#

import netCDF4 as nc
import numpy as np
from scipy.interpolate import griddata
from netCDF4 import Dataset as dat
import matplotlib.pyplot as plt

from funcs_tide import *


grdname = 'grid_rotated_SUL_2_NEST2_smaler.nc';
ininame = 'prooceano_myocean_mais_13_03_ini.nc';

nudname = 'son_2.nc';

G = get_ROMS_grid(grdname);

spherical = G.hgrid.spherical;

[Lr,Mr] = G.hgrid.lat_rho.shape;

Nr = len(G.vgrid.s_rho);

# Set switches for state variables to nudge.

LnudgeM2CLM    = True;           # nudging 2D momentum
LnudgeM3CLM    = True;           # nudging 3D momentum
LnudgeTCLM     = True;           #nudging tracers (usually T-S)
LnudgeTgeneric = False;           # nudging generic tracers

# Set NetCDF variables to process.  Initialize inverse nudging
# coefficients with zeros.  Recall that division by zero is not
# defined and we will give "Inf".  Therefore, we just need to set
# only the values in the desired areas and the rest can be zero
# (for no nudging) because the nudging in ROMS is:
#
#      F(...,new) = F(...,new) +
#                   dt * F_nudgcoef * (Fclm - F(...,new))


nudcoef = 1./5.; # Para fazer o nudging em todo o dominio por igual.
VarNUD = [];
F = [];

if spherical:
  spherical=True
else:
  spherical=False

if (spherical):
  VarNUD = ['lon_rho', 'lat_rho']
else:
  VarNUD = ['x_rho', 'y_rho']


if (LnudgeM2CLM):
  VarNUD.append('M2_NudgeCoef');
  M2_NudgeCoef = np.ones([Lr,Mr]) * nudcoef;                # RHO-points

if (LnudgeM3CLM):
  VarNUD.append('M3_NudgeCoef');
  M3_NudgeCoef = np.ones([Nr,Lr,Mr]) * nudcoef;             # RHO-points

if (LnudgeTCLM):
  VarNUD.append('temp_NudgeCoef')
  VarNUD.append('salt_NudgeCoef')
  temp_NudgeCoef = np.ones([Nr,Lr,Mr]) * nudcoef;
  salt_NudgeCoef = np.ones([Nr,Lr,Mr]) * nudcoef;


if (LnudgeTgeneric):
  VarNUD.append('tracer_NudgeCoef');
  tracer_NudgeCoef = zeros([Nr,Lr,Mr]);
  
  
#--------------------------------------------------------------------------
# Create Nudging coefficients NetCDF file: build creation parameters
# structure, S.
#--------------------------------------------------------------------------

S            = [];
Tindex       = [];
ReplaceValue = np.nan;
PreserveType = True;
Unlimited    = False;                   # time dimension is umlimited
nctype       = 'nc_double';             # input data is in double precision

ncfile = dat(nudname, mode='w', clobber='true', format='NETCDF3_CLASSIC')
ncfile.createDimension('lon', size=Mr)
ncfile.createDimension('lat', size=Lr)
ncfile.createDimension('s_rho', size=Nr)

ncfile.createVariable('spherical', 'i')
setattr(ncfile.variables['spherical'], 'long_name', 'Grid type logical switch')
setattr(ncfile.variables['spherical'], 'flag_values', [0,1])
setattr(ncfile.variables['spherical'], 'flag_meanings', 'Cartesian spherical')
ncfile.variables['spherical'][:]  = int(spherical)

ncfile.createVariable('lon_rho', 'd', dimensions=('lat', 'lon'))
setattr(ncfile.variables['lon_rho'], 'long_name', 'longitude of RHO-points')
setattr(ncfile.variables['lon_rho'], 'units', 'degree_east')
setattr(ncfile.variables['lon_rho'], 'standard_name', 'longitude')
ncfile.variables['lon_rho'][:]  =  G.hgrid.lon_rho

ncfile.createVariable('lat_rho', 'd', dimensions=('lat', 'lon'))
setattr(ncfile.variables['lat_rho'], 'long_name', 'latitude of RHO-points')
setattr(ncfile.variables['lat_rho'], 'units', 'degree_north')
setattr(ncfile.variables['lat_rho'], 'standard_name', 'latitude')
ncfile.variables['lat_rho'][:]  =  G.hgrid.lat_rho

ncfile.createVariable(VarNUD[2], 'd', dimensions=('lat', 'lon'))
setattr(ncfile.variables[VarNUD[2]], 'long_name', '2D momentum inverse nudging coefficients')
setattr(ncfile.variables[VarNUD[2]], 'units', 'day-1')
#setattr(ncfile.variables[VarNUD[2]], 'standard_name', 'latitude')
ncfile.variables[VarNUD[2]][:]  =  M2_NudgeCoef

ncfile.createVariable(VarNUD[3], 'd', dimensions=('s_rho', 'lat', 'lon'))
setattr(ncfile.variables[VarNUD[3]], 'long_name', '3D momentum inverse nudging coefficients')
setattr(ncfile.variables[VarNUD[3]], 'units', 'day-1')
#setattr(ncfile.variables[VarNUD[3]], 'standard_name', 'latitude')
ncfile.variables[VarNUD[3]][:]  =  M3_NudgeCoef

ncfile.createVariable(VarNUD[4], 'd', dimensions=('s_rho', 'lat', 'lon'))
setattr(ncfile.variables[VarNUD[4]], 'long_name', 'temp inverse nudging coefficients')
setattr(ncfile.variables[VarNUD[4]], 'units', 'day-1')
#setattr(ncfile.variables[VarNUD[4]], 'standard_name', 'latitude')
ncfile.variables[VarNUD[4]][:]  =  temp_NudgeCoef

ncfile.createVariable(VarNUD[5], 'd', dimensions=('s_rho', 'lat', 'lon'))
setattr(ncfile.variables[VarNUD[5]], 'long_name', 'salt inverse nudging coefficients')
setattr(ncfile.variables[VarNUD[5]], 'units', 'day-1')
#setattr(ncfile.variables[VarNUD[5]], 'standard_name', 'latitude')
ncfile.variables[VarNUD[5]][:]  =  salt_NudgeCoef


#inneror = 1./30.;                        # x days at interior limit
inneror = 0;
outeror = 2;                         #  x days at boundary very strong
#outeror = 1./30.;                         #  x days at boundary weak
widthor = 10.;                           #  x points

work  = np.zeros([Nr,Lr,Mr])


#--------------------------------------------------------------------------
# Set inverse time scales.
#--------------------------------------------------------------------------

IstrR = 0;
IendR = int(Lr);
JstrR = 0;
JendR = int(Mr);

allnud=15                                 #number of bottom layer to apply complete nudge
innerall = 0.;                        # x days at interior limit
outerall = 2.;                         #  x days at boundary weak
widthall = 10.;                           #  x points

bottom3dnudge=True                  #activate stronger nudging in the bottom

for k in range(Nr):
  if k in range(allnud):
  
    outer=outerall
    inner=innerall
    width=widthall
 
    work[k,:] = work[k,:] + inner
    
    for i in range(IendR):                     # Eastern Boundary
      for j in range(JendR-int(width),JendR):
        work[k,i,j] = outer + (JendR -1 - j) * (inner - outer) / width;
     
    for i in range(IendR):                     # western Boundary
      for j in range(JstrR,JstrR+int(width)):
        work[k,i,j] = inner + (width - j) * (outer - inner) / width;
      
    h=0      
    for i in range(IendR-int(width),IendR):                     # nothern Boundary
      for j in range(JstrR+int(width)-h,JendR-int(width)+h):
        work[k,i,j] = outer + (IendR -1 - i) * (inner - outer) / width;
      h=h+1

    h=int(width)
    for i in range(IstrR,IstrR+int(width)):                     # southern Boundary
      for j in range(JstrR+int(width)-h,JendR-int(width)+h):
        work[k,i,j] = inner + (width - i) * (outer - inner) / width;
      h=h-1

   
  else:
  
    outer=outeror
    inner=inneror
    width=widthor  
    
    
    work[k,:] = work[k,:] + inner    
    
    for i in range(IendR):                     # Eastern Boundary
      for j in range(JendR-int(width),JendR):
        work[k,i,j] = outer + (JendR -1 - j) * (inner - outer) / width;
     
    for i in range(IendR):                     # western Boundary
      for j in range(JstrR,JstrR+int(width)):
        work[k,i,j] = inner + (width - j) * (outer - inner) / width;
      
    h=0      
    for i in range(IendR-int(width),IendR):                     # nothern Boundary
      for j in range(JstrR+int(width)-h,JendR-int(width)+h):
        work[k,i,j] = outer + (IendR -1 - i) * (inner - outer) / width;
      h=h+1

    h=int(width)
    for i in range(IstrR,IstrR+int(width)):                     # southern Boundary
      for j in range(JstrR+int(width)-h,JendR-int(width)+h):
        work[k,i,j] = inner + (width - i) * (outer - inner) / width;
      h=h-1

if bottom3dnudge:
  vwindow=10
  for i in range(int(Lr)):
    for j in range(int(Mr)):
      work[allnud-vwindow:allnud,i,j]=np.linspace(work[allnud-vwindow,i,j],work[allnud,i,j],vwindow)
    

ncfile.variables[VarNUD[2]][:]  =  work[-1,:]


for i in range(3,len(VarNUD)):
  ncfile.variables[VarNUD[i]][:]  =  work
  VarNUD[i]

  
ncfile.sync()
ncfile.close()
print "Done."
print ""  
  