
#########   SET PATHS !!!! ########################

maindir=`pwd`                           #main directory of the project


getmercator=${maindir}/forecast_myocean_data/get_myocean_best_v2.sh   # path to the file to download mercator

mercatordata=${maindir}/mercator_data/                 # Path to the directory to store Mercator data (tip: create a new)

forocean=${maindir}/rotinas_boundary_initial/mercator/   #Path to the directory with .py rotines to create ocean conditions



getgfs=${maindir}/forecast_gfs_data/get_gfs_manual_netcef.sh   # path to the file to download gfs data

gfsdata=${maindir}/gfs_data/                           # Path to the directory to store GF data (tip: create a new)

foratm=${maindir}/rotinas_boundary_initial/gfs/          #Path to the directory with .py rotines to create meteo conditions

##############################################

numdays=1                          #numdays to run ROMS
numdayshind=1                          #numdays to run ROMS hindcast

dly=0          # If 1 begin simulation 1 day before to accomodate Data assimilation (MUR release data 1 day after)



newini=ocean_in_forecast_ciclone.in                 #name of greatest .in file

hisinterdad=4

hisinterson=4   #hours interval for nested grid (important for nesting)

#son_grids=(azul_son_case1_newBAT_2_Mcu.nc azul_son_case2_newlm.nc)  #name of son grids seperated by space ex: son_grids=(son1.nc son2.nc)

DoNest=TRUE         # TRUE OR FALSE  no space between =, if TRUE perform nesting

#son_grids=(CRONOS_ES_36v2ALI1_grd.nc)  #name of son grids seperated by space ex: 
son_grids=(gz_small_9_7.nc)
#son_grids=(abc2.nc)  #name of son grids seperated by space ex: son_grids=(son1.nc son2.nc)

#DTSON=(150 50)  #timestep in seconds of son grids seperated by space 

DTSON=80

fathergrid='gz_great_9_7.nc'
#fathergrid='CRONOS_ES_36v2ALI1_grd.nc'

#DT=100.          # timestep in second of the father grid
DT=200.          # timestep in second of the father grid

visc2_1=100

visc2_2=100

################################ INFO NEST 1

#hindcastnest=hiscast_nest1.in

rotatedid=True
cuttedid=True
cuttedidM=True

lim=2



theta_b=2.0
theta_s=5.0
tcline=250.
klevels=20
Vtransform=2
Vstretching=4


theta_b_1=2.0
theta_s_1=5.0
tcline_1=250.
klevels_1=20
Vtransform_1=2
Vstretching_1=4

tide_st='2021-05-25 12:00:00'  #the ini time when building tide file

############## NUDGING FILES

NUDGECLIM=FALSE                          #FALSE OR TRUE

dadnud=dad.nc

son_nud_1=son_1.nc

son_nud_2=son_2.nc


######## Coordinates Mercator data

xmin=-61
xmax=-20
ymin=-45
ymax=-10
