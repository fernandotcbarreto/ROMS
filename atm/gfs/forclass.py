from parameters import *
from datetime import date
from netCDF4 import Dataset as dat
import numpy as np


class inputs(object):
   def __init__(self, fname,var_name,L,M,nt,lat,lon,time, index, vec, timename):
      self.fname=fname
      self.L=L
      self.M=M
      self.nt=nt
      self.LNf=lon
      self.LAf=lat
      self.var=var_name
      self.vec=vec
      self.index=index
      self.timevar=time    
      self.timename = timename
      
   def change_var(self,fname,variable, data, timename):
      self.fname=fname
      self.var=variable
      self.vec=data
      self.timename = timename

class ncgene(inputs):      
   def __init__(self, fname,var_name,L,M,nt,lat,lon,time, index, vec,timename):
      super(ncgene,self).__init__(fname,var_name,L,M,nt,lat,lon,time, index, vec, timename)
      
   def start(self):
     ncfile = dat(self.fname, mode='w', clobber='true', format='NETCDF3_CLASSIC')
     
     ncfile.createDimension('lon', size=self.L)
     ncfile.createDimension('lat', size=self.M)
     ncfile.createDimension(self.timename, size=self.nt+1)


     ncfile.createVariable('spherical', 'i')
     setattr(ncfile.variables['spherical'], 'long_name', 'Grid type logical switch')
     setattr(ncfile.variables['spherical'], 'flag_values', [0,1])
     setattr(ncfile.variables['spherical'], 'flag_meanings', 'Cartesian spherical')
     ncfile.variables['spherical'][:]  = 1


     ncfile.createVariable('lon', 'd', dimensions=('lat', 'lon'))
     setattr(ncfile.variables['lon'], 'long_name', 'longitude')
     setattr(ncfile.variables['lon'], 'units', 'degree_east')
     setattr(ncfile.variables['lon'], 'standard_name', 'longitude')
     ncfile.variables['lon'][:]  = self.LNf

     ncfile.createVariable('lat', 'd', dimensions=('lat', 'lon'))
     setattr(ncfile.variables['lat'], 'long_name', 'latitude')
     setattr(ncfile.variables['lat'], 'units', 'degree_north')
     setattr(ncfile.variables['lat'], 'standard_name', 'latitude')
     ncfile.variables['lat'][:]  = self.LAf

################################################################################## TIME VARIABLE

     if self.var == 'Uwind' or self.var == 'Vwind':
       ncfile.createVariable(self.timename, 'd', dimensions=self.timename)
       setattr(ncfile.variables[self.timename], 'long_name', 'surface wind time')
#       setattr(ncfile.variables[self.timename], 'units', 'modified Julian day')
       setattr(ncfile.variables[self.timename], 'units', 'seconds since 2013-01-01 00:00:00')
       setattr(ncfile.variables[self.timename], 'calendar', 'gregorian') # 
       setattr(ncfile.variables[self.timename], 'field', str(self.timename)+', scalar, series')
       setattr(ncfile.variables[self.timename], 'time', self.timename)
       
     elif self.var == 'rain': 
       ncfile.createVariable(self.timename, 'd', dimensions=self.timename)
       setattr(ncfile.variables[self.timename], 'long_name', 'rain fall rate time')
#       setattr(ncfile.variables[self.timename], 'units', 'modified Julian day')
       setattr(ncfile.variables[self.timename], 'units', 'seconds since 2013-01-01 00:00:00')
       setattr(ncfile.variables[self.timename], 'calendar', 'gregorian')
       setattr(ncfile.variables[self.timename], 'field', str(self.timename)+', scalar, series')
       setattr(ncfile.variables[self.timename], 'time', self.timename)
       
     elif self.var == 'lwrad_down': 
       ncfile.createVariable(self.timename, 'd', dimensions=self.timename)
       setattr(ncfile.variables[self.timename], 'long_name', 'solar longwave radiation time')
#       setattr(ncfile.variables[self.timename], 'units', 'modified Julian day')
       setattr(ncfile.variables[self.timename], 'units', 'seconds since 2013-01-01 00:00:00')
       setattr(ncfile.variables[self.timename], 'calendar', 'gregorian')
       setattr(ncfile.variables[self.timename], 'field', str(self.timename)+', scalar, series')
       setattr(ncfile.variables[self.timename], 'time', self.timename)

     elif self.var == 'Tair': 
       ncfile.createVariable(self.timename, 'd', dimensions=self.timename)
       setattr(ncfile.variables[self.timename], 'long_name', 'surface air temperature time')
#       setattr(ncfile.variables[self.timename], 'units', 'modified Julian day')
       setattr(ncfile.variables[self.timename], 'units', 'seconds since 2013-01-01 00:00:00')
       setattr(ncfile.variables[self.timename], 'calendar', 'gregorian')
       setattr(ncfile.variables[self.timename], 'field', str(self.timename)+', scalar, series')
       setattr(ncfile.variables[self.timename], 'time', self.timename)
       
     elif self.var == 'Pair': 
       ncfile.createVariable(self.timename, 'd', dimensions=self.timename)
       setattr(ncfile.variables[self.timename], 'long_name', 'surface air pressure time')
#       setattr(ncfile.variables[self.timename], 'units', 'modified Julian day')
       setattr(ncfile.variables[self.timename], 'units', 'seconds since 2013-01-01 00:00:00')
       setattr(ncfile.variables[self.timename], 'calendar', 'gregorian')
       setattr(ncfile.variables[self.timename], 'field', str(self.timename)+', scalar, series')
       setattr(ncfile.variables[self.timename], 'time', self.timename)
       
       
     elif self.var == 'swrad': 
       ncfile.createVariable(self.timename, 'd', dimensions=self.timename)
       setattr(ncfile.variables[self.timename], 'long_name', 'solar shortwave radiation time')
#       setattr(ncfile.variables[self.timename], 'units', 'modified Julian day')
       setattr(ncfile.variables[self.timename], 'units', 'seconds since 2013-01-01 00:00:00')
       setattr(ncfile.variables[self.timename], 'calendar', 'gregorian')
       setattr(ncfile.variables[self.timename], 'field', str(self.timename)+', scalar, series')
       setattr(ncfile.variables[self.timename], 'time', self.timename)


     elif self.var == 'Qair': 
       ncfile.createVariable(self.timename, 'd', dimensions=self.timename)
       setattr(ncfile.variables[self.timename], 'long_name', 'surface relative humidity time')
#       setattr(ncfile.variables[self.timename], 'units', 'modified Julian day')
       setattr(ncfile.variables[self.timename], 'units', 'seconds since 2013-01-01 00:00:00')
       setattr(ncfile.variables[self.timename], 'calendar', 'gregorian')
       setattr(ncfile.variables[self.timename], 'field', str(self.timename)+', scalar, series')
       setattr(ncfile.variables[self.timename], 'time', self.timename)
###################################################################################################### FORCING VARIABLE

     if self.var == 'Uwind':
       ncfile.createVariable(self.var, 'd', dimensions=(self.timename, 'lat', 'lon'))
       setattr(ncfile.variables[self.var], 'long_name', '10m u-wind component')
       setattr(ncfile.variables[self.var], 'long_name', 'u-wind')
       setattr(ncfile.variables[self.var], 'units', 'meter second-1')
       setattr(ncfile.variables[self.var], 'units', 'm/s')
       setattr(ncfile.variables[self.var], 'field', 'uwind, scalar, series')
       setattr(ncfile.variables[self.var], 'time', self.timename)
       setattr(ncfile.variables[self.var], 'coordinates', 'lon lat '+str(self.timename))
       setattr(ncfile.variables[self.var], 'coordinates', 'lon lat '+str(self.timename))
 
     elif self.var == 'Vwind':
       ncfile.createVariable(self.var, 'd', dimensions=(self.timename, 'lat', 'lon'))
       setattr(ncfile.variables[self.var], 'long_name', '10m v-wind component')
       setattr(ncfile.variables[self.var], 'long_name', 'v-wind')
       setattr(ncfile.variables[self.var], 'units', 'meter second-1')
       setattr(ncfile.variables[self.var], 'units', 'm/s')
       setattr(ncfile.variables[self.var], 'field', 'vwind, scalar, series')
       setattr(ncfile.variables[self.var], 'time', self.timename)
       setattr(ncfile.variables[self.var], 'coordinates', 'lon lat '+str(self.timename))
       setattr(ncfile.variables[self.var], 'coordinates', 'lon lat '+str(self.timename))
       
     elif self.var == 'rain':
       ncfile.createVariable(self.var, 'd', dimensions=(self.timename, 'lat', 'lon'))
       setattr(ncfile.variables[self.var], 'long_name', 'rain fall rate')
       setattr(ncfile.variables[self.var], 'long_name', 'rain fall rate')
       setattr(ncfile.variables[self.var], 'units', 'kilogram meter-2 second-1')
       setattr(ncfile.variables[self.var], 'units', 'kilogram meter-2 second-1')
       setattr(ncfile.variables[self.var], 'field', 'rain, scalar, series')
       setattr(ncfile.variables[self.var], 'time', self.timename)
       setattr(ncfile.variables[self.var], 'coordinates', 'lon lat '+str(self.timename))
       setattr(ncfile.variables[self.var], 'coordinates', 'lon lat '+str(self.timename))
       
     elif self.var == 'lwrad_down':
       ncfile.createVariable(self.var, 'd', dimensions=(self.timename, 'lat', 'lon'))
       setattr(ncfile.variables[self.var], 'long_name', 'downward longwave radiation flux')
       setattr(ncfile.variables[self.var], 'long_name', 'downward longwave radiation flux')
       setattr(ncfile.variables[self.var], 'units', 'Watts meter-2')
       setattr(ncfile.variables[self.var], 'units', 'Watts meter-2')
       setattr(ncfile.variables[self.var], 'positive', 'downward flux, heating')
       setattr(ncfile.variables[self.var], 'positive', 'downward flux, heating')
       setattr(ncfile.variables[self.var], 'negative', 'upward flux, cooling')
       setattr(ncfile.variables[self.var], 'negative', 'upward flux, cooling')
       setattr(ncfile.variables[self.var], 'field', 'lwrad, scalar, series')
       setattr(ncfile.variables[self.var], 'time', self.timename)
       setattr(ncfile.variables[self.var], 'coordinates', 'lon lat '+str(self.timename))
       setattr(ncfile.variables[self.var], 'coordinates', 'lon lat '+str(self.timename))

     elif self.var == 'Tair':
       ncfile.createVariable(self.var, 'd', dimensions=(self.timename, 'lat', 'lon'))
       setattr(ncfile.variables[self.var], 'long_name', 'surface air temperature')
       setattr(ncfile.variables[self.var], 'long_name', 'surface air temperature')
       setattr(ncfile.variables[self.var], 'units', 'Celsius')
       setattr(ncfile.variables[self.var], 'units', 'Celsius')
       setattr(ncfile.variables[self.var], 'field', 'Tair, scalar, series')
       setattr(ncfile.variables[self.var], 'time', self.timename)
       setattr(ncfile.variables[self.var], 'coordinates', 'lon lat '+str(self.timename))
       setattr(ncfile.variables[self.var], 'coordinates', 'lon lat '+str(self.timename))
       
     elif self.var == 'Pair':
       ncfile.createVariable(self.var, 'd', dimensions=(self.timename, 'lat', 'lon'))
       setattr(ncfile.variables[self.var], 'long_name', 'surface air pressure')
       setattr(ncfile.variables[self.var], 'long_name', 'surface air pressure')
       setattr(ncfile.variables[self.var], 'units', 'milibar')
       setattr(ncfile.variables[self.var], 'units', 'milibar')
       setattr(ncfile.variables[self.var], 'field', 'Pair, scalar, series')
       setattr(ncfile.variables[self.var], 'time', self.timename)
       setattr(ncfile.variables[self.var], 'coordinates', 'lon lat '+str(self.timename))
       setattr(ncfile.variables[self.var], 'coordinates', 'lon lat '+str(self.timename))
       
       
     elif self.var == 'swrad':
       ncfile.createVariable(self.var, 'd', dimensions=(self.timename, 'lat', 'lon'))
       setattr(ncfile.variables[self.var], 'long_name', 'solar shortwave radiation')
       setattr(ncfile.variables[self.var], 'long_name', 'shortwave radiation')
       setattr(ncfile.variables[self.var], 'units', 'Watts meter-2')
       setattr(ncfile.variables[self.var], 'units', 'Watts meter-2')
       setattr(ncfile.variables[self.var], 'positive', 'downward flux, heating')
       setattr(ncfile.variables[self.var], 'positive', 'downward flux, heating')
       setattr(ncfile.variables[self.var], 'negative', 'upward flux, cooling')
       setattr(ncfile.variables[self.var], 'negative', 'upward flux, cooling')
       setattr(ncfile.variables[self.var], 'field', 'swrad, scalar, series')
       setattr(ncfile.variables[self.var], 'time', self.timename)
       setattr(ncfile.variables[self.var], 'coordinates', 'lon lat '+str(self.timename))
       setattr(ncfile.variables[self.var], 'coordinates', 'lon lat '+str(self.timename))
       
     elif self.var == 'Qair':
       ncfile.createVariable(self.var, 'd', dimensions=(self.timename, 'lat', 'lon'))
       setattr(ncfile.variables[self.var], 'long_name', 'relative humidity')
       setattr(ncfile.variables[self.var], 'long_name', 'relative humidity')
       setattr(ncfile.variables[self.var], 'units', 'percentage')
       setattr(ncfile.variables[self.var], 'units', 'percentage')
       setattr(ncfile.variables[self.var], 'field', 'Qair, scalar, series')
       setattr(ncfile.variables[self.var], 'time', self.timename)
       setattr(ncfile.variables[self.var], 'coordinates', 'lon lat '+str(self.timename))
       setattr(ncfile.variables[self.var], 'coordinates', 'lon lat '+str(self.timename))
       

     setattr(ncfile, 'title', ROMS_title)
     setattr(ncfile, 'grd_file', grdname)

     i_dtref=self.index;

     self.vec=np.transpose(self.vec,(2,1,0))

     ind=0

     for t in range(i_dtref[0][0],len(self.timevar)):
       ui=np.squeeze(self.vec[t,::])
       ncfile.variables[self.timename][ind]  = self.timevar[t]
       ncfile.variables[self.var][ind,:,:]  = ui
       ind=ind+1
      
     ncfile.sync()
     ncfile.close()
     print ("Done.")
     print ("")