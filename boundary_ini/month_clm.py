#monthly  climatology

import matplotlib.dates as dates
import numpy as np
from netCDF4 import Dataset
from parameters_bry_in import *



yearclm=2013

monthiclm=1                        #initial month

monthfclm=5                        #final month, will be computded 3 months


##################### parameters


datei=str(yearclm)+'-'+str(monthiclm).zfill(2)+'-'+str(1).zfill(2) + ' ' +str(0).zfill(2)+':'+str(0).zfill(2)+':'+str(0).zfill(2)

datef=str(yearclm)+'-'+str(monthfclm).zfill(2)+'-'+str(1).zfill(2) + ' ' +str(0).zfill(2)+':'+str(0).zfill(2)+':'+str(0).zfill(2)

a=[datei, datef]

numdays= np.diff(dates.datestr2num(a)) #goes until last day of monthfclm - 1

n_time=int(numdays[0])

a = dates.datestr2num(a)

time_vec = np.arange(numdays) + a[0]

n_time = np.arange(monthfclm - monthiclm) + 1

nc_hycom=input_path + dates.num2date(time_vec[0]).strftime("%Y%m%d")+'.nc'

infofile=Dataset(nc_hycom)
info=infofile['uo'][:]
uclm = np.zeros( (len(n_time),) + info.shape[1:])
vclm = np.zeros( (len(n_time),) + info.shape[1:])
tempclm = np.zeros( (len(n_time),) + info.shape[1:])
saltclm = np.zeros( (len(n_time),) + info.shape[1:])
sshclm =  np.zeros( (len(n_time),) + info.shape[2:])
numdaysclm=np.zeros(len(n_time),)
bry_time = np.zeros(len(n_time),)
check_max = np.zeros(len(n_time),)


numdaysclm = numdaysclm*0
for clml in range(len(n_time)): 
  #load HYCOM DATA 
  for rt in range(numdays):
    if (dates.num2date( time_vec[rt]).month == n_time[clml]):
      nc_hycom=input_path + dates.num2date( time_vec[rt]).strftime("%Y%m%d")+'.nc'
      numdaysclm[clml] = numdaysclm[clml] + 1
      file=Dataset(nc_hycom)
      uclm[clml] = uclm[clml] + np.squeeze(file['uo'][:]).filled()
      vclm[clml] = vclm[clml] + np.squeeze(file['vo'][:]).filled()
      maxvalue = np.squeeze(file['vo'][:]).filled().max()
      if maxvalue > check_max[clml]:
        check_max[clml] = maxvalue
      tempclm[clml] = tempclm[clml] + np.squeeze(file['thetao'][:]).filled()
      saltclm[clml] = saltclm[clml] + np.squeeze(file['so'][:]).filled()
      sshclm[clml] = sshclm[clml]+ np.squeeze((file['zos'][:])).filled()
  uclm[clml] = uclm[clml]/numdaysclm[clml] 
  vclm[clml] = vclm[clml]/numdaysclm[clml]
  tempclm[clml] = tempclm[clml]/numdaysclm[clml]
  saltclm[clml] = saltclm[clml]/numdaysclm[clml]
  sshclm[clml] = sshclm[clml]/numdaysclm[clml]
  
uclm[uclm<-100] = np.nan
vclm[vclm<-100] = np.nan
tempclm[tempclm<-100] = np.nan
saltclm[saltclm<-100] = np.nan
sshclm[sshclm<-100] = np.nan


bry_time = numdaysclm/2.
for i in range(len(n_time)-1):
  datei = str(yearclm)+'-'+str(int(n_time[i])).zfill(2)+'-'+str(int(numdaysclm[i])).zfill(2) + ' ' +str(0).zfill(2)+':'+str(0).zfill(2)+':'+str(0).zfill(2)
  bry_time[i+1] = bry_time[i+1] + np.diff(dates.datestr2num([timeref,datei])) 
 
bry_time[0]  = bry_time[0] - 1  #since the timecount starts at day 1

bry_time = bry_time * 24*60*60

numdays=len(n_time)


      print(nc_hycom)
      print(n_time[clml])
  