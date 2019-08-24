function TH4_Bat, filename,liste=liste, level= level, Lambert = Lambert, quiet=quiet, save=save, $
s_channel=s_channel,b_channel=b_channel, direct=direct, powel=powel, delta=delta

;+
; NAME
;   TH4_Bat
;
; PURPOSE
;  Fit reflectance and T from radiance observations of 67P by Virtis H and M
;     Fit under continuity constraint (NO), + *NO* multithread version (?)
;   Use MPfitfun
;
; CALLING SEQUENCE:
;   therm = th4_bat(!dirpreland+'MTP009/STP025/CALIBRATED/I1_00372813212.CAL',level=3,/s_,/lambert,/quiet,/direct ,/save)
;
; INPUTS:
;   filename: provides path and name of individual file in dir (.CAL) just like for virtispds
; OR
;   nothing if the keyword /liste is activated
;
; OUTPUTS:
;   a structure with the output variables 
;   
; KEYWORDS:
;     Quiet : minimizes print
;     liste : allow to select several files in a directory (only M OR H) [it would be easily possible to adapt the programm to change VH or VM case at each iteration]
;     level : choose the level (level 1 : first derivation of the temperature, 2 : T+radiance factor, 3 : T of first step, radiance factor of second step and last derivation of T [which does not work well for the moment])
;     lambert : if activated, Lambert law will be used, otherwise : Lommel Seeliger. /!\ LS is not available at step 2 because it is needed to have the more complicated expression of the law by solving an integral
;     save : save the structure for each processed cube
;     s_channel : use the channel of Stephane (in case of VIRTIS-H)
;     b_channel : use the channel of Batiste (in case of VIRTIS-H)
;     direct : direct computation of the radiance factor at the second step (no use of MPFITFUN which provide same result taking much more time)
;     powel : use another technique to do the fit of the radiance at the first step, provide approx. the same result as MPFITFUN
;     delta : use another technique at the third step (minimization of the chi² in place of the radiance)
;     
; COMMENTS:
;  Calls invtherm_mh4_bat in a loop
;  Need to use/select the cube .CAL => radiance. For H AND M. For M, the artefact-removed version of the cube "*.QUB_rad_ar_1.0.0.dat" will be open and used for the process
;
; EXAMPLES:
;
; HISTORY:
;   From Tlutetia :
;     SE, LESIA, March 2009: Tsteins2
;     SE, LESIA, Oct 2009: use more recent GEO file + Now uses ENVI reflectance (instead of radiance)
;          + this version used only a selection of channels to fit
;     Tlutetia:
;     SE, LESIA, Nov 2011: adapted for Lutetia, with more generic function v_invtherm
;                         Interactive tuning is detailed in Therm_invert.jou
;   Tlutetia2 :
;     SE, LESIA, Jan 2013: use function v_invtherm3, fit with continuity constraint
;     First try with a short loop to test batch mode on Lesia04 (and Nadja)
;     v_invtherm3 will now return upon error with rf =(-1) (cst vector)
;     Writes intermediate results in dummy files (to be erased manually)
;   Tlutetia2_split :
;     SE, LESIA, Jan 2013: from Tlutetia_split and Tlutetia2
;   Tlutetia4_split :
;     SE, LESIA, Feb 2013: from Tlutetia2_split
;   Tlutetia4c_split :
;     SE, LESIA, Feb 2013: from Tlutetia4_split
;   Tlutetia4d_split :
;     SE, LESIA, 21 Feb 2013: from Tlutetia4c_split
;   TH3_split :
;     SE, LESIA, 5 Oct 2014: from Tlutetia4d_split + v_mapstp
;   TH3 :
;     SE, LESIA, 12 Jan 2016: from TH3_split, simplified (no split loop) + Fixed correspondance of data and angles
;     5 Feb 2016: quick fix of rf when no reflectance returned
;   Modification by B. Rousseau - IPAG - 06/2018
;     - programme adapté pour VIRTIS-M
;     - divers arrangements pour changer l'entrée, les sorties, etc 

if keyword_set(quiet) then quiet=1 else quiet=0
Lambert = keyword_set(Lambert) 
if (get_login_info()).(0) eq 'DELLOU' then savepath = 'V:\2018_06_Thermique\Out\TH4_'
if (get_login_info()).(0) eq 'IPAGOU' then savepath = 'C:\Users\rousseab\Desktop\Taff_local\2018_06_Thermique\Out2\TH4_'

VH=(VM=0)
if ~keyword_set(liste) then toprocess = filename else toprocess = dialog_pickfile(path=!dirdata,/multiple,filter=['T1*.CAL','I1*.CAL']) ; if liste is activated we want to process several cube, this line allow to chose them
if strmatch(toprocess(0),'*I1_*') eq 1 then VM=1
if strmatch(toprocess(0),'*T1_*') eq 1  then VH=1

Nf = N_elements(toprocess)

; ouverture des lambda selon VM ou VH
if VM eq 1 then begin
  lam = v_lamm(/ir)/1000.
  ; lecture du spectre solaire
  container=read_ascii(!virtis_sun+'KURUCZ_VIRTIS_RESAMPLED_IR_HIGH_v10.asc',data_start=7,header=entete)
  RSoleil=double(reform(container.field1(2,*))) ;extract solar irradiance from container
  if ~quiet then print, "Reading solar spectrum: ", !virtis_sun+'KURUCZ_VIRTIS_RESAMPLED_IR_HIGH_v10.asc'
endif
if VH eq 1 then begin
  lam = v_lamh(/ros)
  merge = merge_cube()
  if keyword_set(s_channel) then can = MERGE.S_TOOLS.S_CHANNEL
  if keyword_set(b_channel) then can = MERGE.B_TOOLS.MY_CHANNEL
  if ~keyword_set(b_channel) and ~keyword_set(s_channel) then begin
    print,'Channel keyword is missing, S_channel has been activated'
    can = MERGE.S_TOOLS.S_CHANNEL
  endif
  lam = lam(can)
  ; lecture du spectre solaire
  RSoleil = fltarr(432,8)
  fich = 'HsoleilRos2014.asc'     ; only for H
  Fich = FilePath(fich, root='/Users/fandrieu/Documents/Programmes/Archive_Batiste/'+'Virtis_Soleil\')
  openr, 1, Fich
  readf, 1, RSoleil
  close, 1
  if ~quiet then print, "Reading solar spectrum: ", Fich
endif

; début de l'itération sur les différents cubes
For i=0, Nf-1 do begin
  filename =toprocess[i]
  print, 'Process '+filename
  
  cal = virtispds(filename, /si) ; open the cube
  if VM eq 1 then radm = find_dat(cal,/radiance)
  geocube = find_geo2(cal) ; find the corresponding geocube
  
  ; date and distance to the Sun in AU 
  start_date = v_pdspar(cal.label, 'START_TIME')
  dista = v_pdspar(cal.label, 'SOLAR_DISTANCE') / 149.598E6
  
  ; cosinus incidence and emergence, work for M and H ; resolution
  mu = reform(cos(geocube.qube(11,*,*)/10000. * !dtor)) ; emergence (OK for H and M)
  mu0 = reform(cos(geocube.qube(10,*,*)/10000.* !dtor)) ; incidence (OK for H and M)
  resol = reso(geocube)
  
  ;Look for pixels not on the nucleus (OK for H and M)
  notnucleus = where(reform(geocube.qube(17,*,*)*geocube.qube_coeff(17)) GE 100., Nbpix)
  sz = cal.qube_dim
  
  if VH eq 1 then begin
    spec2 = cal.qube
    nelement = sz(1)
  endif
  if VM eq 1 then begin
    x = sz(1)
    y = sz(2)
    nelement = x*y
    spec2 = reform(radm,432,x*y)
    mu = reform(mu,x*y)
    mu0 = reform(mu0,x*y)
  endif
  spec2(*,notnucleus)='NaN' ; fill the location not in the nucleus by NaN, to avoid process at next step but to keep the same dimension
    
;  if ~quiet then print, "You must compile invtherm_h4_mp before in order to compile the fit function"
  t0 = systime(1)

  For ii=0, nelement-1 do begin
    if ~quiet then print, ii
    
    testnan = where(finite(spec2(*,ii)),/null) ; 'cause it's full of NaN : two possibilities => we are in the sky or this is due the .DAT file of VM (first column often)
    if testnan eq !null then begin
      T1=(T3=(chi2_1=(chi2_3=(niter1=(niter3=(rf1='NaN'))))))
      if level eq 3 then begin
        rf3 = dblarr(n_elements(lam)) 
        rf3(*)='NaN'        
      endif
      goto,next
    endif
    
    invtherm_mh4_bat, spec2(*,ii), mu0(ii), mu(ii), dista, lam,T1, T3,chi2_1,chi2_3,rf1,rf3, Sflux=Rsoleil, $ ; init is not used anymore 
      Niter1= Niter1,Niter3= Niter3, lambert=lambert, level=level, quiet=quiet, s_channel=s_channel,b_channel=b_channel, $
      direct=direct, powel=powel, delta=delta, vm=vm, vh=vh

    next:    
    If N_elements(rf1result) eq 0 then rf1result=rf1 else rf1result = [rf1result,rf1] ; facteur de radiance level 1
    if N_elements(T1result) EQ 0 then T1result = T1 else T1result = [T1result,T1] ; T1 = Température déterminé au level 1
    if N_elements(Chiresult1) EQ 0 then Chiresult1=chi2_1 else Chiresult1=[Chiresult1, chi2_1]
    if N_elements(Nbiter1) EQ 0 then Nbiter1=niter1 else Nbiter1=[Nbiter1, float(niter1)] ; float needed because the NaN are a problem otherwise
    if level eq 3 then begin
      if N_elements(rf3result) eq 0 then rf3result=rf3 else rf3result = [[rf3result],[rf3]] ; facteur de radiance level 1
      if N_elements(T3result) EQ 0 then T3result = T3 else T3result = [T3result,T3] ; T3 = Température déterminé au level 3
      if N_elements(Chiresult3) EQ 0 then Chiresult3=chi2_3 else Chiresult3=[Chiresult3, chi2_3]
      if N_elements(Nbiter3) EQ 0 then Nbiter3=niter3 else Nbiter3=[Nbiter3, float(niter3)]
    endif
  endfor

  Texec=systime(1) -T0
  print, '  Tp exec (s): ', Texec

  if VM then begin
    rf1result = reform(rf1result,x,y)
    T1result = reform(T1result,x,y)
    Chiresult1 = reform(Chiresult1,x,y)
    Nbiter1 = reform(Nbiter1,x,y)
    mu0 = reform(mu0,x,y)
    mu = reform(mu,x,y)
    if level eq 3 then begin
      rf3result = reform(rf3result,432,x,y)
      T3result = reform(T3result,x,y)
      chiresult3 = reform(CHIRESULT3,x,y)
      Nbiter3 = reform(Nbiter3,x,y)
    endif
  endif
  
  MTP_name = strmid(filename,strpos(filename,'MTP'),6)
  STP_name = strmid(filename,strpos(filename,'STP'),6)
  if lambert then corr='Lambert' else corr='LS'
  outfilename = MTP_name+'_'+STP_name+'_'+corr+'_'+strsplit(v_pdspar(cal.label,'PRODUCT_ID'),'"',/extract)
  
  if level ne 3 then begin
    T3result=(chi3=(rf3result=(Nbiter3='NaN')))
  endif  
  
  out = {wavel:lam,$ ; lambda
    cosinci:mu0,$
    cosemer:mu,$
    T1:float(T1result),$ ; température level 1
    T3:float(T3result),$ ; température level 1
    Chi2_1:float(chiresult1),$
    Chi2_3:float(chiresult3),$
    rf1:float(rf1result),$ ; facteur de radiance level 1
    rf3:float(rf3result),$ ; facteur de radiance level 3
    Nbiter1:Nbiter1,$ ; nombre d'itération level 1
    Nbiter3:Nbiter3,$ ; nombre d'itération level 3
    Lambert:Lambert,$
    filename:filename,$
    direct:direct,$
    powel:powel,$
    level:level,$
    start_date:start_date}
  
  T1=(T3=(Chi2_1=(Chi2_3=(rf1=(rf3=(Nbiter1=(Nbiter3=!NULL)))))))
  
   if keyword_set(save) then begin
     outfile = savepath+outfilename+'.sav'
    save,out,filename=outfile
   endif


  ; for the next iteration
  rf1result=(T1result=(Chiresult1=(Nbiter1=(rf3result=(T3result=(Chiresult3=(Nbiter3=!NULL)))))))

endfor

return,out

end