maindir=/home/fernando/roms/src/Projects/operational   #main directory of the project
mercatordata=${maindir}/mercator_data/

forocean='/home/fernando/roms/src/Projects/operational/rotinas_boundary_initial/mercator/'

numdays=7                          #numdays to run ROMS

newini=ocean_in_forecast_ciclone.in                 #name of greatest .in file



hisinterson=6   #hours interval for nested grid (important for nesting)

#son_grids=(azul_son_case1_newBAT_2_Mcu.nc azul_son_case2_newlm.nc)  #name of son grids seperated by space ex: son_grids=(son1.nc son2.nc)

son_grids=(grid_rotated_SUL_2_NEST_smaler.nc grid_rotated_SUL_2_NEST2_smaler.nc)  #name of son grids seperated by space ex: son_grids=(son1.nc son2.nc)


DTSON=(150 150)  #timestep in seconds of son grids seperated by space 

#DTSON=150

fathergrid='grid_rotated_SUL_2_smaler.nc'
################################ INFO NEST 1

#hindcastnest=hiscast_nest1.in

rotatedid=True
cuttedid=True

lim=0


theta_b=0.4
theta_s=5.0
tcline=100.
klevels=40
Vtransform=2
Vstretching=4


############## NUDGING FILES

NUDGECLIM=TRUE

dadnud=dad.nc

son_nud_1=son_1.nc

son_nud_2=son_2.nc


######## Coordinates Mercator data

xmin=-52
xmax=-23
ymin=-32
ymax=-14