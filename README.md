smr2613
=======

The tutorial codes used in "School on Parallel Programming and Parallel Architecture for HPC and Developer School for HPC applications in Earth Sciences + Symposium on HPC and Data-Intensive Applications in Earth Sciences"

The detailed information about the event can be found in the following link:

http://indico.ictp.it/event/a13229/

Installation
============
* The Makefile can be used to install the test application (just type make command). In this case, ESMFMKFILE environment variable must be set to point out esmf.mk file in the ESMF installation directory.
* The test application requires "atm_grd.nc" and "ocn_grd.nc" files to define model grids in the ESMF side. To create those files user need to use "grid2scrip.ncl" NCL script under proc/ directory. In this case, user could change the atm_grid_file and ocn_grid_file parameters (or variables) in the NCL script. Also note that ESMF regriiding capability is supported only NCL versions >= 6.2.0.
* The next step is to edit configuration file (namelist.rc) of the test application. In this case, tile information for each component and grid definitions (created in previous step) must be defined correctly. The total number of processor used in the application is equal to the multipication of number of tiles in each direction (i.e. 2 x 6 = 12).
* Then main.job job submission script can be used to run the test application.
* The application creates bunch of netCDF files for exchanged variables, masking etc. The application basically demostrates two step interpolation: (1) Interpolation that performs regridding only ober ocean using bilinear type interpolation, (2) Extrapolation that used to fill the unmapped grid points close to the coastline (due to the mismatch between land-sea mask of the model componantes and their horizontal resolutions).
* The fix.sh shell script might be used to add coordinate information to the ESMF generated netCDF files using NCO commands.
* The error comes from interpolation and extrapolation can be visualizaed using plot_err.ncl NCL script under proc/ directory. In this case, the exchanged fields are compared with their analytical solution to calculate the error.
