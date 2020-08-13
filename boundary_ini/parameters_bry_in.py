
#number of columns/lines to trim the input file in the interpolation.
#Example: if lim = 2, the code will employ 4 gridlines for each interpolation.
      
lim=2

######################################################################################

noextra=True            # in vertical 1D interpolation uses the nearest value in extrapolation
                        # if False continues with the linear interpolation

timeini = '2020-06-26 00:00:00'

timeend = '2020-07-03 00:00:00'

timeref = '2013-01-01 00:00:00'                   # Reference time at bry file. Azul project is seconds since 2013-01-01 00:00:00]

input_path='/home/fernando/roms/src/Projects/operational/mercator_data/MYOCEAN_AZUL_FORECAST_'

#input_path='C:\Users\Fernando\Downloads\CMEMS_BEST_ANALYSIS_20130213_20130513.tar\CMEMS_BEST_ANALYSIS_20130213_20130513\CMEMS_BEST_ANALYSIS_'    # directory with myocean data


run_name = '20200626'      # name for bry and ini files

fname_grd = 'grid_rotated_ciclone.nc'                       # ROMS grid   

#fname_grd = 'rtte.nc'


## Stretching curve parameters.
theta_b = 0.4
theta_s = 5.0
tcline = 100.
klevels = 40
Vtransform = 2
Vstretching = 4


Spherical = True

#monthly  climatology

yearclm=2016

monthiclm=4                        #initial month

monthfclm=5                        #final month, will be computded 3 months

monthly=False
##################### parameters


########NESTING PARAMETERS

rotated=True     # Activate if FATHER grid is rotated (son is not a problem)
cutted=True       # Ativate if you want to interpolate from strips of father grid (pay attention when father grid is rottated)

fname_grd_son='grid_rotated_ciclone_NEST.nc'   #son grid

coup_files=['HIS_FILE_20200626_5D0-20200703_5D0_hind.nc']

run_name_son='prooceano_myocean_ciclone_NEST'
