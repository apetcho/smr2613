;-----------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;-----------------------------------------------------------
;***********************************************************
;*** Ufuk Turuncoglu ***************************************
;*** Send bug reports to turuncu@be.itu.edu.tr *************
;***********************************************************
begin
  ;--- set constants ---
  missing = 1.0e20
  pi = 4.0*atan(1.0)
  deg2rad = pi/180.0
 
  ;--- open files ---
  nc1 = addfile("remap_1_forward.nc", "r")
  nc2 = addfile("remap_2_forward.nc", "r")

  ;--- read variables ---
  var1 = nc1->data
  var1 = mask(var1, var1 .lt. missing/2.0, True) 
  var2 = nc2->data
  var2 = mask(var2, var2 .lt. missing/2.0, True) 

  lat2d = nc1->lat_rho
  lon2d = nc1->lon_rho

  ;--- calculate analytical solution ---
  ;--- pseudo spherical harmonics (L=32, M=16) ---
  var0 = 2.0+sin(2.0*lat2d*deg2rad)^16*cos(16.0*lon2d*deg2rad)
  var0 = mask(var0, ismissing(var2), False)
 
  ;--- create workspace for plot ---
  wks = gsn_open_wks("newpdf", "plot_err")
  gsn_define_colormap(wks, "MPL_jet")
  p = NhlNewColor(wks, 0.8, 0.8, 0.8)

  ;--- set plot resources ---
  res = True
  res@gsnFrame = False
  res@gsnLeftString = ""
  res@gsnRightString = ""
  res@gsnSpreadColors = True
  res@gsnSpreadColorEnd = -3
  res@cnFillOn = True
  res@cnFillMode = "RasterFill"
  res@cnInfoLabelOn = False
  res@cnLinesOn = False
  res@cnLineLabelsOn = False
  res@lbLabelBarOn = False
  res@pmTickMarkDisplayMode = "Always"
  res@tiXAxisFontHeightF = 0.008
  res@tiYAxisFontHeightF = 0.008
  res@tmXBLabelFontHeightF = 0.008
  res@tmYLLabelFontHeightF = 0.008
  res@gsnAddCyclic = False

  res@tmXBOn = False
  res@tmXBLabelsOn = False
  res@tmXTOn = False
  res@tmXTLabelsOn = False
  res@tmYROn = False
  res@tmYRLabelsOn = False
  res@tmYLOn = False
  res@tmYLLabelsOn = False

  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels = fspan(1.0, 3.0, 21)

  ;--- mapping properties --- 
  res@mpDataBaseVersion = "HighRes"
  res@mpProjection = "LambertConformal"
  res@mpOutlineDrawOrder = "PostDraw"
  res@mpGridAndLimbOn = False
  minlon = min(lon2d)
  maxlon = max(lon2d)
  onesix = (maxlon-minlon)/6.0
  std1 = minlon+onesix
  std2 = maxlon-onesix
  print("[warning] -- one-six-rule standard parallels ("+\
  sprintf("%6.2f", std1)+","+sprintf("%6.2f", std2)+") are used.")
  res@mpLambertParallel1F = std1
  res@mpLambertParallel2F = std2

  res@mpLimitMode = "Corners"
  res@mpLeftCornerLatF = min(lat2d)-0.25
  res@mpLeftCornerLonF = min(lon2d)-0.15
  res@mpRightCornerLatF = max(lat2d)+0.25
  res@mpRightCornerLonF = max(lon2d)+0.50
  clon = (max(lon2d)+min(lon2d))/2.0
  res@mpLambertMeridianF = clon
  print("[warning] -- default center longtitude ("+sprintf("%6.2f", clon)+") is used. ")

  dims = dimsizes(lat2d)
  nlat = dims(0)
  nlon = dims(1)
  fact = tofloat(nlon)/tofloat(nlat)
  delete(dims)
  res@vpHeightF = 0.2
  res@vpWidthF = fact*res@vpHeightF

  res@tmYLOn = True
  res@tmYLLabelsOn = True
  res@tiYAxisString = "Latitude"

  res@vpXF = 0.1
  res@vpYF = 0.9
  var0@lat2d = lat2d
  var0@lon2d = lon2d
  plot1 = gsn_csm_contour_map(wks, var0, res)

  res@tmYLOn = False
  res@tmYLLabelsOn = False
  res@tiYAxisString = ""

  res@lbLabelBarOn = True
  res@lbLabelFontHeightF = 0.008
  res@lbOrientation = "Vertical"
  res@lbTitleOffsetF = 0.001
  res@pmLabelBarWidthF = 0.06

  res@vpXF = 0.1+res@vpWidthF*0.98
  res@vpYF = 0.9
  var1@lat2d = lat2d
  var1@lon2d = lon2d
  plot2 = gsn_csm_contour_map(wks, var1, res)

  res@tmXBOn = True
  res@tmXBLabelsOn = True
  res@tmYLOn = True
  res@tmYLLabelsOn = True
  res@tiYAxisString = "Latitude"
  res@tiXAxisString = "Longitude"
  res@lbLabelBarOn = False

  res@vpXF = 0.1
  res@vpYF = 0.9-res@vpHeightF
  var2@lat2d = lat2d
  var2@lon2d = lon2d
  plot3 = gsn_csm_contour_map(wks, var2, res)

  gsn_define_colormap(wks, "MPL_seismic")
  p = NhlNewColor(wks, 0.8, 0.8, 0.8)
  delete(res@cnLevels)
  res@cnLevels = fspan(-0.25, 0.25, 21)

  res@tmYLOn = False
  res@tmYLLabelsOn = False
  res@tiYAxisString = ""

  res@lbLabelBarOn = True
  res@lbLabelFontHeightF = 0.008
  res@lbOrientation = "Vertical"
  res@lbTitleOffsetF = 0.001
  res@pmLabelBarWidthF = 0.06

  res@vpXF = 0.1+res@vpWidthF*0.98
  res@vpYF = 0.9-res@vpHeightF
  ;diff = var2-var0
  diff = (var2/var0)-1.0
  print("min = "+sprintf("%10.5f", min(diff)))
  print("avg = "+sprintf("%10.5f", avg(diff)))
  print("max = "+sprintf("%10.5f", max(diff)))
  diff@lat2d = lat2d
  diff@lon2d = lon2d
  plot4 = gsn_csm_contour_map(wks, diff, res)

  ;--- add title ---
  txres = True
  txres@gsnFrame = False
  txres@txJust = "CenterLeft"
  txres@txPerimOn = True
  txres@txFontHeightF = 0.012
  txres@txBackgroundFillColor = "white"

  title = (/ "FUNC", "REMAP_1", "REMAP_2", "(REMAP_2/FUNC)-1" /)
  do j = 0, 3
    txres@txPosXF = 0.1+tofloat(j%2)*res@vpWidthF*0.98+0.015
    txres@txPosYF = 0.9-tofloat(j/2)*res@vpHeightF-0.02
    txid = gsn_create_text(wks, title(j), txres)
    draw(txid)
  end do

  frame(wks)
end
