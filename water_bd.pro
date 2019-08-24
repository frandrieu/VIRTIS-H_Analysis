;;+
; NAME:
;     water_BD
;
; PURPOSE:
;     Compute the water band depth at 3.077µm for a given VIRTIS-H sectrum.
;     Classical linear interpolation between 2.55µm and 3.9µm (biggest water band bounds in VIRTIS-H range)
;
; CALLING SEQUENCE:
;     BD_water=water_bd( spectrum_VH )
;
; INPUTS:
;     spectrum_VH:  any VIRTIS-H specrum
;
; OUTPUTS:
;     BD_water:     water band depth at 3.077µm
;
; KEYWORD PARAMETERS:
;
; EXAMPLE:
;
; LIMITATIONS:
;
; MODIFICATION HISTORY:
;   F ANDRIEU, LESIA, 2017-11-02: Creation
;      Version 1.0


function water_BD, spectrum_VH 

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


  Rc=median(spectrum_VH[pwc])
  Rl=median(spectrum_VH[pwl])
  Rr=median(spectrum_VH[pwr])
  b=( band_center-band_left ) / ( band_right-band_left)
  a=1.-b

  BD_water= 1. - Rc / ( a*Rl + b*Rr )
  return, BD_water


end