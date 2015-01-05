;-----------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
;-----------------------------------------------------------
begin
  ;--- parameters ---
  atm_grid_file = "MED50_DOMAIN000.nc"
  ;atm_grid_file = "TR10_DOMAIN000.nc"
  ocn_grid_file = "BSEA_grd_v1.nc"

  ;--- open atm model file ---
  atm_nc = addfile(atm_grid_file, "r")
  atm_lat2d = atm_nc->xlat
  atm_lon2d = atm_nc->xlon
  atm_msk2d = atm_nc->mask

  ;--- open ocn model grid file ---
  ocn_nc = addfile(ocn_grid_file, "r")
  ocn_lat2d = ocn_nc->lat_rho
  ocn_lon2d = ocn_nc->lon_rho
  ocn_msk2d = ocn_nc->mask_rho

  ;--- create SCRIP definition of grids ---
  ofile = "atm_grd.nc"
  if (.not. isfilepresent(ofile)) then
    opt = True
    opt@ForceOverwrite = True
    opt@PrintTimings = True
    opt@Mask2D = where(atm_msk2d .lt. 2, 1, 0)
    opt@Title = "RegCM grid"
    curvilinear_to_SCRIP(ofile, atm_lat2d, atm_lon2d, opt)
    delete(opt)
  end if

  ;--- create SCRIP definition of destination grid ---
  ofile = "ocn_grd.nc"
  if (.not. isfilepresent(ofile)) then
    opt = True
    opt@ForceOverwrite = True
    opt@PrintTimings = True
    opt@Mask2D = where(ocn_msk2d .gt. 0, 1, 0)
    opt@Title = "ROMS grid"
    curvilinear_to_SCRIP(ofile, ocn_lat2d, ocn_lon2d, opt)
    delete(opt)
  end if
end
