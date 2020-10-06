#modify netCDF4
import numpy as np

from netCDF4 import Dataset 
import matplotlib.pyplot as plt
from roms_setup import run_setup, rho2uvp, get_metrics, spheric_dist
from roms_setup import get_angle, add_topo, process_mask, uvp_mask, smoothgrid
from scipy.interpolate import griddata, interp1d
from roms_setup import get_angle, add_topo, process_mask, uvp_mask, smoothgrid
import netCDF4
from bathy_smoother.LP_bathy_smoothing import LP_smoothing_rx0



gridin='tete.nc'
file=Dataset(gridin)
print('Converting grid' + gridin)

lat_rho=np.real(np.real(file['lat_rho'][:]))

lon_rho=np.real(np.real(file['lon_rho'][:]))

h=np.real(np.real(file['h'][:]))

hraw=np.real(np.real(file['hraw'][:]))

#angle=np.real(np.real(file['angle'][:]))

f=np.real(np.real(file['f'][:]))

pm=np.real(np.real(file['pm'][:]))

pn=np.real(np.real(file['pn'][:]))

dndx=np.real(np.real(file['dndx'][:]))

dmde=np.real(np.real(file['dmde'][:]))

mask_rho=np.real(np.real(file['mask_rho'][:]))

xl=file['xl'][:]

el=file['el'][:]


Lonu, Lonv, Lonp  = rho2uvp(lon_rho)
Latu, Latv, Latp  = rho2uvp(lat_rho)
[masku, maskv, maskp] = uvp_mask(mask_rho) 

angle = get_angle(Latu, Lonu)
#pm, pn, dndx, dmde = get_metrics(Latu, Lonu, Latv, Lonv)

#f0 = 4 * np.pi * np.sin( np.pi * lat_rho/180. ) / ( 24.*3600. )

ETOPO1=True
if ETOPO1:
  print('bathymetry from ETOPO1')
  etopo1=netCDF4.Dataset('etopo1_bed_g2.nc')
  latbat=np.ma.filled(etopo1['lat'][:])
  lonbat=np.ma.filled(etopo1['lon'][:])
  batgt=np.ma.filled(etopo1['topo'][:]).astype(float)
  lonbat,latbat=np.meshgrid(lonbat,latbat)
  batgt=batgt*-1
#  batgt[batgt<0]=np.nan
#  batgt=np.ma.masked_invalid(batgt)
  row,col=np.where((latbat<-18) & (latbat>-40) & (lonbat>-50) & (lonbat<-35))
  row=np.unique(row)
  col=np.unique(col)
  col,row=np.meshgrid(col,row)
  latbat=latbat[row,col]
  lonbat=lonbat[row,col]
  batgt=batgt[row,col]
#  plt.pcolor(lonbat,latbat, batgt, cmap=cmocean.cm.cmap_d['deep']);plt.show() 
  h=griddata((lonbat.ravel(),latbat.ravel()),batgt.ravel(),(lon_rho.ravel(),lat_rho.ravel())).reshape(lat_rho.shape)
#################

ETOPO2=False
if ETOPO2:
  print('bathymetry from ETOPO2')
  h = add_topo(lon_rho, lat_rho, pm, pn,'ETOPO2v2g_f4.nc')   #etopo

#plt.pcolor(lon_rho, lat_rho, h);plt.show()


dad_interp=True

if dad_interp:
  p_file=netCDF4.Dataset('CF_tmz_g026.nc')
  h_parent=p_file['h'][:]
  lon_parent=p_file['lon_rho'][:]
  lat_parent=p_file['lat_rho'][:]
  h_parent=griddata((lon_parent.ravel(),lat_parent.ravel()),h_parent.ravel(),(lon_rho.ravel(),lat_rho.ravel())).reshape(lon_rho.shape)
  nans=np.isnan(h_parent)
  h_parent[nans]=griddata((lon_rho[~nans],lat_rho[~nans]),h_parent[~nans],(lon_rho[nans],lat_rho[nans]), 'nearest')


maskr = h*0
maskr=np.abs(maskr)
maskr[ np.where(h > 10) ] = 1 
maskr = process_mask(maskr)
[masku, maskv, maskp] = uvp_mask(maskr) 

hmin=10

print(angle)


AmpConst=np.zeros(maskr.shape)+10000

SignConst=np.zeros(maskr.shape)
#signindex=np.where((hraw>0)&(hraw<1000))
#SignConst[signindex]=-1

hraw=h.copy()

h=LP_smoothing_rx0(maskr,hraw,0.26,SignConst,AmpConst)

#plt.plot(-h[-1,:]);plt.show()


h = smoothgrid(h, maskr, hmin,hmin,0.26, 1,1)                               # OLD FILTER


hh = h.copy()

hh=h_parent.copy()                        ###MESMA BATI DA MALHA PAI
hraw=h_parent.copy()                        ###MESMA BATI DA MALHA PAI

#############################################################  weights for transition between the two smoothed grids
[Lr,Mr] = lon_rho.shape;

inner = 1;
outer = 0.;                         
width = 20.;                          

work  = np.zeros([Lr,Mr]) + inner;


#--------------------------------------------------------------------------
# Set inverse time scales.
#--------------------------------------------------------------------------

IstrR = 0;
IendR = int(Lr);
JstrR = 0;
JendR = int(Mr);

for i in range(IendR):                     # Eastern Boundary
    for j in range(JendR-int(width),JendR):
      work[i,j] = outer + (JendR -1 - j) * (inner - outer) / width;

for i in range(IendR):                     # western Boundary
    for j in range(JstrR,JstrR+int(width)):
      work[i,j] = inner + (width - j) * (outer - inner) / width;

h=0      
for i in range(IendR-int(width),IendR):                     # nothern Boundary
    for j in range(JstrR+int(width)-h,JendR-int(width)+h):
      work[i,j] = outer + (IendR -1 - i) * (inner - outer) / width;
    h=h+1

h=int(width)
for i in range(IstrR,IstrR+int(width)):                     # southern Boundary
    for j in range(JstrR+int(width)-h,JendR-int(width)+h):
      work[i,j] = inner + (width - i) * (outer - inner) / width;
    h=h-1

##########################################################################################

#plt.plot(-hh[-1,:]);plt.show()
#plt.plot(-h_parent[-1,:]);plt.show()

hson=hh.copy()
h=hh.copy()
#h_parent=hh.copy()          #testando primeiro com tudo uniforme

Lonr=lon_rho.copy()

for i in range(Lonr.shape[0]):
  for j in range(Lonr.shape[1]):
    h[i,j]=(work[i,j]*hson[i,j]) + ((1. - work[i,j])*h_parent[i,j])

#plt.plot(-h[-1,:]);plt.show()



maskr = h*0
maskr=np.abs(maskr)
maskr[ np.where(h > 20.) ] = 1 
maskr = process_mask(maskr)
[masku, maskv, maskp] = uvp_mask(maskr) 

mask_rho =maskr.copy()


plt.pcolor(lon_rho, lat_rho, mask_rho);plt.colorbar();
plt.contour(lon_rho,lat_rho,h,levels=[200],colors=('red'))

plt.show()

plt.pcolor(lon_rho, lat_rho, h);plt.colorbar();plt.show()


plt.pcolor(Lonv, Latv, maskv);plt.colorbar();plt.show()


M, L = Latp.shape


print (' \n' + '==> ' + '  WRITING NETCDF GRID FILE ...\n' + ' ')
import datetime as dt

now = dt.datetime.now()
Lp = L + 1
Mp = M + 1

#if run.spherical == 1:
spherical = 'T'
#else:
#	spherical = 'F'

ncfile = Dataset('dudu_2.nc', mode='w',
    clobber='true', format='NETCDF3_CLASSIC')

# creating DIMENSIONS
ncfile.createDimension('xi_psi', size=L)
ncfile.createDimension('xi_rho', size=Lp)
ncfile.createDimension('xi_u', size=L)
ncfile.createDimension('xi_v', size=Lp)
ncfile.createDimension('eta_psi', size=M)
ncfile.createDimension('eta_rho', size=Mp)
ncfile.createDimension('eta_u', size=Mp)
ncfile.createDimension('eta_v', size=M)
#ncfile.createDimension('bath', size=None)
ncfile.createDimension('one', size=1)
ncfile.createDimension('two', size=2)
ncfile.createDimension('four', size=4)


# creating GLOBAL ATTRIBUTES
setattr(ncfile, 'type', 'new')
setattr(ncfile, 'title', 'grids')
setattr(ncfile, 'history', str(now))


# creating VARIABLES, ATTRIBUTES and ASSIGNING VALUES

# ---------------------------------------------------------------------------
ncfile.createVariable('spherical', 'c')
setattr(ncfile.variables['spherical'], 'long_name', 'Grid type logical switch')
setattr(ncfile.variables['spherical'], 'option_T', 'spherical')
setattr(ncfile.variables['spherical'], 'option_F', 'cartesian')
ncfile.variables['spherical'][:]  = spherical

# ---------------------------------------------------------------------------
ncfile.createVariable('xl', 'd', dimensions=('one'))
setattr(ncfile.variables['xl'], 'long_name', 'domain length in XI-direction')
setattr(ncfile.variables['xl'], 'units', 'meter')
ncfile.variables['xl'][:]  = xl

# ---------------------------------------------------------------------------
ncfile.createVariable('el', 'd', dimensions=('one'))
setattr(ncfile.variables['el'], 'long_name', 'domain length in ETA-direction')
setattr(ncfile.variables['el'], 'units', 'meter')
ncfile.variables['el'][:]  = el

# ---------------------------------------------------------------------------
ncfile.createVariable('hraw', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['hraw'], 'long_name', 'Working bathymetry at RHO-points')
setattr(ncfile.variables['hraw'], 'units', 'meter')
setattr(ncfile.variables['hraw'], 'coordinates', 'lon_rho lat_rho bath')
ncfile.variables['hraw'][:]  = hraw

# ---------------------------------------------------------------------------
ncfile.createVariable('h', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['h'], 'long_name', 'Final bathymetry at RHO-points')
setattr(ncfile.variables['h'], 'units', 'meter')
setattr(ncfile.variables['h'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['h'][:]  = h

# ---------------------------------------------------------------------------
ncfile.createVariable('f', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['f'], 'long_name', 'Coriolis parameter at RHO-points')
setattr(ncfile.variables['f'], 'units', 'second-1')
setattr(ncfile.variables['f'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['f'][:]  = f

# ---------------------------------------------------------------------------
ncfile.createVariable('pm', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['pm'], 'long_name', 'Curvilinear coordinate metric in XI')
setattr(ncfile.variables['pm'], 'units', 'meter-1')
setattr(ncfile.variables['pm'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['pm'][:]  = pm

# ---------------------------------------------------------------------------
ncfile.createVariable('pn', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['pn'], 'long_name', 'Curvilinear coordinate metric in ETA')
setattr(ncfile.variables['pn'], 'units', 'meter-1')
setattr(ncfile.variables['pn'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['pn'][:]  = pn

# ---------------------------------------------------------------------------
ncfile.createVariable('dndx', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['dndx'], 'long_name', 
	'XI derivative of inverse metric factor pn')
setattr(ncfile.variables['dndx'], 'units', 'meter')
setattr(ncfile.variables['dndx'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['dndx'][:]  = dndx

# ---------------------------------------------------------------------------
ncfile.createVariable('dmde', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['dmde'], 'long_name', 
	'ETA derivative of inverse metric factor pm')
setattr(ncfile.variables['dmde'], 'units', 'meter')
setattr(ncfile.variables['dmde'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['dmde'][:]  = dmde

# ---------------------------------------------------------------------------
ncfile.createVariable('lon_rho', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['lon_rho'], 'long_name', 'longitude of RHO-points')
setattr(ncfile.variables['lon_rho'], 'units', 'degree east')
ncfile.variables['lon_rho'][:]  = lon_rho

# ---------------------------------------------------------------------------
ncfile.createVariable('lat_rho', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['lat_rho'], 'long_name', 'latitude of RHO-points')
setattr(ncfile.variables['lat_rho'], 'units', 'degree north')
ncfile.variables['lat_rho'][:]  = lat_rho

# ---------------------------------------------------------------------------
ncfile.createVariable('lon_psi', 'd', dimensions=('eta_psi', 'xi_psi'))
setattr(ncfile.variables['lon_psi'], 'long_name', 'longitude of PSI-points')
setattr(ncfile.variables['lon_psi'], 'units', 'degree east')
ncfile.variables['lon_psi'][:]  = Lonp

# ---------------------------------------------------------------------------
ncfile.createVariable('lat_psi', 'd', dimensions=('eta_psi', 'xi_psi'))
setattr(ncfile.variables['lat_psi'], 'long_name', 'latitude of PSI-points')
setattr(ncfile.variables['lat_psi'], 'units', 'degree north')
ncfile.variables['lat_psi'][:]  = Latp

# ---------------------------------------------------------------------------
ncfile.createVariable('lon_u', 'd', dimensions=('eta_u', 'xi_u'))
setattr(ncfile.variables['lon_u'], 'long_name', 'longitude of U-points')
setattr(ncfile.variables['lon_u'], 'units', 'degree east')
ncfile.variables['lon_u'][:]  = Lonu

# ---------------------------------------------------------------------------
ncfile.createVariable('lat_u', 'd', dimensions=('eta_u', 'xi_u'))
setattr(ncfile.variables['lat_u'], 'long_name', 'latitude of U-points')
setattr(ncfile.variables['lat_u'], 'units', 'degree north')
ncfile.variables['lat_u'][:]  = Latu

# ---------------------------------------------------------------------------
ncfile.createVariable('lon_v', 'd', dimensions=('eta_v', 'xi_v'))
setattr(ncfile.variables['lon_v'], 'long_name', 'longitude of V-points')
setattr(ncfile.variables['lon_v'], 'units', 'degree east')
ncfile.variables['lon_v'][:]  = Lonv

# ---------------------------------------------------------------------------
ncfile.createVariable('lat_v', 'd', dimensions=('eta_v', 'xi_v'))
setattr(ncfile.variables['lat_v'], 'long_name', 'latitude of V-points')
setattr(ncfile.variables['lat_v'], 'units', 'degree north')
ncfile.variables['lat_v'][:]  = Latv

# ---------------------------------------------------------------------------
ncfile.createVariable('angle', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['angle'], 'long_name', 'angle between XI-axis and EAST')
setattr(ncfile.variables['angle'], 'units', 'radians')
setattr(ncfile.variables['angle'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['angle'][:]  = angle

# ---------------------------------------------------------------------------
ncfile.createVariable('mask_rho', 'd', dimensions=('eta_rho', 'xi_rho'))
setattr(ncfile.variables['mask_rho'], 'long_name', 'mask on RHO-points')
setattr(ncfile.variables['mask_rho'], 'flag_values', '0, 1')
setattr(ncfile.variables['mask_rho'], 'flag_meanings', 'land, water')
setattr(ncfile.variables['mask_rho'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['mask_rho'][:]  = mask_rho

# ---------------------------------------------------------------------------
ncfile.createVariable('mask_u', 'd', dimensions=('eta_u', 'xi_u'))
setattr(ncfile.variables['mask_u'], 'long_name', 'mask on U-points')
setattr(ncfile.variables['mask_u'], 'flag_values', '0, 1')
setattr(ncfile.variables['mask_u'], 'flag_meanings', 'land, water')
setattr(ncfile.variables['mask_u'], 'coordinates', 'lon_u lat_u')
ncfile.variables['mask_u'][:]  = masku

# ---------------------------------------------------------------------------
ncfile.createVariable('mask_v', 'd', dimensions=('eta_v', 'xi_v'))
setattr(ncfile.variables['mask_v'], 'long_name', 'mask on V-points')
setattr(ncfile.variables['mask_v'], 'flag_values', '0, 1')
setattr(ncfile.variables['mask_v'], 'flag_meanings', 'land, water')
setattr(ncfile.variables['mask_v'], 'coordinates', 'lon_v lat_v')
ncfile.variables['mask_v'][:]  = maskv

# ---------------------------------------------------------------------------
ncfile.createVariable('mask_psi', 'd', dimensions=('eta_psi', 'xi_psi'))
setattr(ncfile.variables['mask_psi'], 'long_name', 'mask on PSI-points')
setattr(ncfile.variables['mask_psi'], 'flag_values', '0, 1')
setattr(ncfile.variables['mask_psi'], 'flag_meanings', 'land, water')
setattr(ncfile.variables['mask_psi'], 'coordinates', 'lon_rho lat_rho')
ncfile.variables['mask_psi'][:]  = maskp

ncfile.sync()

print (' \n' + '==> ' + '  ############################################  ...\n' + ' ')
print (' \n' + '==> ' + '        GRID FILE SUCCESSFULLY CREATED          ...\n' + ' ')
print (' \n' + '==> ' + '  ############################################  ...\n' + ' ')




























