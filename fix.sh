#!/bin/bash

ncks -v lon_rho,lat_rho,mask_rho proc/BSEA_grd_v1.nc grid.nc
ncks -a -A grid.nc remap_1_forward.nc
ncks -a -A grid.nc remap_2_forward.nc
ncks -a -A grid.nc gcomp2_data_exp.nc
ncks -a -A grid.nc gcomp2_data_imp.nc
rm -f grid.nc

ncks -v xlon,xlat,mask proc/MED50_DOMAIN000.nc grid.nc
#ncks -v xlon,xlat,mask proc/TR10_DOMAIN000.nc grid.nc
ncks -a -A grid.nc remap_1_backward.nc 
ncks -a -A grid.nc remap_2_backward.nc 
ncks -a -A grid.nc gcomp1_data_exp.nc
ncks -a -A grid.nc gcomp1_data_imp.nc
rm -f grid.nc
