;;+
; NAME:
;     WATER_DETECT_VH_ALL_CUBES
;
; PURPOSE:
;     Adaptation of WATER_BD_VH_CUBE to look for water in the whole VIRTIS-H database
;     Classical linear interpolation between 2.55µm and 3.9µm (biggest water band bounds in VIRTIS-H range)
;
; CALLING SEQUENCE:
;    WATER_DETECT_VH_ALL_CUBES, LIST_CUBES, DETECT_CANDIDATES
;
; INPUTS:
;     LIST_CUBES:  A list of VIRTIS-H CALIBRATED cubes
;
; OUTPUTS:
;     DETECT_CANDIDATES:  a structure containing 3
;
; KEYWORD PARAMETERS:
;     BASELINE=BASELINE: ads a baseline of the specified value to the data to avoid division by 0
;     CRITERION        : the value for the band depth above which a stectrum is kept for further 
;                       inverstigation. If no value is specified, it is set to 0.5.
;                
;     
; EXAMPLE:
;
; LIMITATIONS:
;
; MODIFICATION HISTORY:
;   F ANDRIEU, LESIA, 2017-11-02: Creation
;      Version 1.0


pro WATER_DETECT_VH_ALL_CUBES, LIST_CUBES, DETECT_CANDIDATES, BASELINE=BASELINE, CRITERION=CRITERION

;    restore, '/Users/fandrieu/Documents/Save_Macbook_pro/Documents/Inversion_Achille-dir/Const_optiques/index_ice_rf_warren_2008.txt-ascii_template'
;    data=read_ascii('/Users/fandrieu/Documents/Save_Macbook_pro/Documents/Inversion_Achille-dir/Const_optiques/index_ice_rf_warren_2008.txt', template=temp)
;    lambda=data.FIELD1
;    n=data.FIELD2
;    k=data.FIELD3
;    w=where((lambda gt 1.5) and (lambda lt 5.5))
;    lambda_select=lambda[w]
;    n_select=n[w]
;    k_select=k[w]

  restore,'/Users/fandrieu/Documents/Programmes/IDL_sav_files/VIRTIS-H_waves.sav'
  restore,'/Users/fandrieu/Documents/Programmes/IDL_sav_files/VIRTIS-H_orders.sav'
  ;band_center=lambda_select[where(k_select eq max(k_select))]
  ;water_band_center=band_center[0]
  ;pw=where(abs(wavelengths-water_band_center
  ;) lt 0.05 and orders eq 3, /null)

  

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
  if (keyword_set(baseline)) then base=baseline else base=0.
  if (keyword_set(criterion)) then crit=criterion else crit=0.5
  
  strong_bd=!NULL
  names=!NULL
  pos=!NULL
  n_spectrum_detect=0
  n_cubes_detect=0
  for cube_index=0, n_elements(list_cubes)-1 do begin
      
    cube = virtispds(list_cubes[cube_index], /silent)
    spectra=cube.qube
    n_sp=cube.QUBE_DIM[1]
    Rc=fltarr(n_sp)
    Rl=fltarr(n_sp)
    Rr=fltarr(n_sp)
    BD_water=fltarr(n_sp)
    
    Rc=median(spectra[pwc,*], dimension=1)+base
    Rl=median(spectra[pwl,*], dimension=1)+base
    Rr=median(spectra[pwr,*], dimension=1)+base
  
    BD_water= 1. - Rc / ( a*Rl + b*Rr )
    sb=where(BD_water gt crit, /NULL)
    if (n_elements(sb) ne 0) then begin
      strong_bd=[strong_BD, BD_water[sb]]
      names=[names, replicate(list_cubes[cube_index],n_elements(sb))]
      pos=[pos, sb]
      n_cubes_detect=n_cubes_detect+1
      n_spectrum_detect=n_spectrum_detect+n_elements(sb)
    endif
    
    
  endfor
  print,n_spectrum_detect, '  candidates found in ', n_cubes_detect, ' VIRTIS-H cubes'
  if (n_spectrum_detect eq 0) then DETECT_CANDIDATES=!NULL else $
    DETECT_CANDIDATES={CUBE_NAME:names, ACQUISITION_NUMBER:pos, BAND_DEPTH:strong_bd}
  return


end