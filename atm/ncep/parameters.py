######################################################################
#
# romstools_param: common parameter file for the preprocessing
#                  of ROMS simulations using ROMSTOOLS
#
#                  This file is used by make_grid.m, make_forcing.m, 
#                  make_clim.m, make_biol.m, make_bry.m, make_tides.m,
#                  make_NCEP.m, make_OGCM.m, make_...
# 
#  Further Information:  
#  http://www.brest.ird.fr/Roms_tools/
#  
#  This file is part of ROMSTOOLS
#
#  ROMSTOOLS is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published
#  by the Free Software Foundation; either version 2 of the License,
#  or (at your option) any later version.
#
#  ROMSTOOLS is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#  MA  02111-1307  USA
#
#  Copyright (c) 2005-2006 by Patrick Marchesiello and Pierrick Penven 
#  e-mail:Pierrick.Penven@ird.fr  
#
#  Updated    6-Sep-2006 by Pierrick Penven
#  Updated    2006/10/05 by Pierrick Penven (add tidegauge observations)
#  Updated    24-Oct-2006 by Pierrick Penven (diagnostics, chla etc...)
#
######################################################################
#
# 1 General parameters
#
######################################################################
#
#  ROMS title names and directories
#
#ROMS_title  = 'STS';
#ROMS_config = 'STS';
#ROMSTOOLS_dir = 'C:\Users\Public\Documents\matlab\Roms_tools';
#
# ROMS file names (grid, forcing, bulk, climatology, initial)
#
#grdname=['Y:\PROJETOS\04_P&D\ROMS\3D\MODELO_AZUL\AZUL_112\grade\azul_grd2_orig.nc'];

# Objective analysis decorrelation scale [m]
# (if Roa=0: simple extrapolation method; crude but much less costly)
#
#Roa=300e3;
#Roa=0;
#
#interp_method = 'cubic';           # Interpolation method: 'linear' or 'cubic'
#
#makeplot     = 0;                 # 1: create a few graphics after each preprocessing step
#
######################################################################
#
# 2 Grid parameters
#   used by make_grid.m (and others..)
#
#

######################################################################
#
# 6 Temporal parameters (used for make_tides, make_NCEP, make_OGCM)
#
######################################################################
#
Yorig         = 2013;               # reference time for reference time Roms 
                               


input_path='C:/Users/Fernando/Desktop/prooceano_roms/atm_ncep_reanalysis_2/2013_01_12/'

output_path = 'C:/Users/Fernando/Desktop/prooceano_roms/atm_ncep_reanalysis_2/output/'