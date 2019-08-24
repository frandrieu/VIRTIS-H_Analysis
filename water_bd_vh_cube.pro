;;+
; NAME:
;     WATER_BD_VH_CUBE
;
; PURPOSE:
;     Compute the water band depth at 3.077µm for every sectrum in a given VIRTIS-H cube.
;     Classical linear interpolation between 2.55µm and 3.9µm (biggest water band bounds in VIRTIS-H range)
;
; CALLING SEQUENCE:
;     WATER_BD_VH_CUBE, cube_VH, BD_water
;
; INPUTS:
;     CUBE_VH:  any VIRTIS-H CALIBRATED CUBE
;
; OUTPUTS:
;     BD_water:     water band depth at 3.077µm for every spectrum in the cube, stored in an array.
;
; KEYWORD PARAMETERS:
;     BASELINE=BASELINE: ads a baseline of the specified value to the data to avoid division by 0
; EXAMPLE:
;
; LIMITATIONS:
;
; MODIFICATION HISTORY:
;   F ANDRIEU, LESIA, 2017-11-02: Creation
;      Version 1.0


pro WATER_BD_VH_CUBE, CUBE_VH, BD_WATER, BASELINE=BASELINE

  ;  restore, '/Users/fandrieu/Documents/Save_Macbook_pro/Documents/Inversion_Achille-dir/Const_optiques/index_ice_rf_warren_2008.txt-ascii_template'
  ;  data=read_ascii('/Users/fandrieu/Documents/Save_Macbook_pro/Documents/Inversion_Achille-dir/Const_optiques/index_ice_rf_warren_2008.txt', template=temp)
  ;  lambda=data.FIELD1
  ;  n=data.FIELD2
  ;  k=data.FIELD3
  ;  w=where((lambda gt 1.5) and (lambda lt 5.5))
  ;  lambda_select=lambda[w]
  ;  n_select=n[w]
  ;  k_select=k[w]
 
  restore,'/Users/fandrieu/Documents/Programmes/IDL_sav_files/VIRTIS-H_waves.sav'
  restore,'/Users/fandrieu/Documents/Programmes/IDL_sav_files/VIRTIS-H_orders.sav'
  ;band_center=lambda_select[where(k_select eq max(k_select))]
  ;water_band_center=band_center[0]
  ;pw=where(abs(wavelengths-water_band_center) lt 0.05 and orders eq 3, /null)

  cube = virtispds(cube_VH, /silent)
  
  spectra=cube.qube
  band_center=3.077
  order_center=3
  band_left=2.55
  order_left=4
  band_right=3.9
  order_right=1
  water_band_center=band_center[0]
  pwc=where(abs(wavelengths-band_center) lt 0.05 and orders eq order_center, /null)
  pwl=where(abs(wavelengths-band_left) lt 0.05 and orders eq order_left, /null)
  pwr=where(abs(wavelengths-band_right) lt 0.05 and orders eq order_right, /null)
  b=( band_center-band_left ) / ( band_right-band_left)
  a=1.-b

  n_sp=cube.QUBE_DIM[1]
  Rc=fltarr(n_sp)
  Rl=fltarr(n_sp)
  Rr=fltarr(n_sp)
  BD_water=fltarr(n_sp)
  if (keyword_set(baseline)) then base=baseline else base=0.
  Rc=median(spectra[pwc,*], dimension=1)+base
  Rl=median(spectra[pwl,*], dimension=1)+base
  Rr=median(spectra[pwr,*], dimension=1)+base
  
  BD_water= 1. - Rc / ( a*Rl + b*Rr )
  return


end