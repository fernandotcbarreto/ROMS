
#########   SET PATHS !!!! ########################

maindir=/home/fernando/roms/src/Projects/operational   #main directory of the project


getmercator='/home/fernando/roms/src/Projects/forecast_myocean_data/get_myocean_best_v2.sh'   # path to the file to download mercator

mercatordata=${maindir}/mercator_data/                 # Path to the directory to store Mercator data (tip: create a new)

forocean='/home/fernando/roms/src/Projects/operational/rotinas_boundary_initial/mercator/'   #Path to the directory with .py rotines to create ocean conditions



getgfs='/home/fernando/roms/src/Projects/forecast_gfs_data/get_gfs_manual_netcef.sh'   # path to the file to download gfs data

gfsdata=${maindir}/gfs_data/                           # Path to the directory to store GF data (tip: create a new)

foratm='/home/fernando/roms/src/Projects/operational/rotinas_boundary_initial/gfs/'          #Path to the directory with .py rotines to create meteo conditions

##############################################

numdays=7                          #numdays to run ROMS

newini=ocean_in_forecast_ciclone.in                 #name of greatest .in file


hisinterdad=12

hisinterson=6   #hours interval for nested grid (important for nesting)

#son_grids=(azul_son_case1_newBAT_2_Mcu.nc azul_son_case2_newlm.nc)  #name of son grids seperated by space ex: son_grids=(son1.nc son2.nc)

DoNest=TRUE         # TRUE OR FALSE  no space between =, if TRUE perform nesting

son_grids=(abc1.nc abc2.nc)  #name of son grids seperated by space ex: son_grids=(son1.nc son2.nc)


DTSON=(300 150)  #timestep in seconds of son grids seperated by space 

#DTSON=150

fathergrid='azul_grd2.nc'

DT=360.          # timestep in second of the father grid
################################ INFO NEST 1

#hindcastnest=hiscast_nest1.in

rotatedid=True
cuttedid=True
cuttedidM=False

lim=2


theta_b=0.4
theta_s=5.0
tcline=100.
klevels=40
Vtransform=2
Vstretching=4


############## NUDGING FILES

NUDGECLIM=FALSE                          #FALSE OR TRUE

dadnud=dad.nc

son_nud_1=son_1.nc

son_nud_2=son_2.nc


######## Coordinates Mercator data

xmin=-52
xmax=-23
ymin=-32
ymax=-14