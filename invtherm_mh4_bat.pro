Pro invtherm_mh4_bat, input_radiance, m0, m, dista,lam,T1, T3,chi2_1,chi2_3,rf1,rf3, temp = temp, refmoy= refmoy, Sflux=Sflux, $
  Niter1= Niter1,Niter3= Niter3, Lambert = Lambert, level = level, Quiet= quiet, s_channel=s_channel,b_channel=b_channel, $
  direct=direct, powel=powel, delta=delta, vm=vm, vh=vh
  common param1,  uncert, deltaR, Lcrack
  common param2, flux1, lams
  common param7, Radiance
  common res, fun
  common param5, RSoleil, distS, mu, mu0
  common param6, TempF
  
; NAME
;   invtherm_mh4_bat
;
; PURPOSE
;  Fit reflectance and T from radiance observations by Virtis H and VIRTIS M
;     Uses Virtis measured radiance
; 	  First tests for H on 67P, uses merged spectra (no order overlap)
;	first fit T and average refl, then spectral reflectance ("radiance factor") only, then T again (in progress...work but not work in fact)
;
; CALLING SEQUENCE:
;   called by th4_bat
;
; INPUTS:
;   input_radiance: measured radiance (for H the .CAL cube, for M the derived-data .DAT cube) => see TH4_BAT
;   mu0, mu : cosines of angles of incidence and emergence
;   dista : Sun distance in AU
;   lam: wavelengths vector, /!\ to the ascending or descending order depending of VM or VH. The program is adapt to taking into account this 
;   
; OUTPUTS:
;   T1 : temperature output at the first step
;   T3 : temperature output at the third step
;   chi2_1 and chi2_3 : CHI² for 1st and 3rd steps but do they make sense ? not sure... 
;   rf1 and rf3 : radiance factor (it's called spectral reflectance sometimes...) from 1st and 3rd steps 
;
; KEYWORDS:
;     temp : initial temperature (default is 200K) ; not used
;     refmoy : initial reflectance (default is to estimate from shorter wavelengths) ; not used
;     Lambert : Use Lambert function instead of the default Lommell-Seeliger 
;     Quiet : minimizes print and plot
;     Level: steps performed. 1: fit T only; 2: Reflectance from fit T; 3: T from fit reflectance (iterated)
;     Sflux: Solar flux
;     Niter1= Niter1,
;     Niter3= Niter3, 
;     s_channel : use the channel of Stephane (in case of VIRTIS-H)
;     b_channel : use the channel of Batiste (in case of VIRTIS-H)
;     direct : direct computation of the radiance factor at the second step (no use of MPFITFUN which provide same result taking much more time)
;     powel : use another technique to do the fit of the radiance at the first step, provide approx. the same result as MPFITFUN
;     delta : use another technique at the third step (minimization of the chi² in place of the radiance)
;     vm, vh : precise VM or VH is in use
;
; COMMENTS:
;  Uses solar spectrum computed from Kurucz 1997, only spectrum covering this range
;  In case of pb check mu & mu0 (should be cosines!) +  Sun distance +lamref + v_lamH (set for Rosetta)
;  Lommel-Seeliger fit available in option (change fit routine in call to fitting routine)
;  
;  Note by B.Rousseau : 
;   - the step 3 is now working but does not provide better fit, and the temperature is the same as derived in the first step. Then, T is probably overestimated in most case 
;       The problem coulb be due to the second step which provide the radiance factor by fitting the radiance curve, but no matter the way to do this (direct computation or fit with MPFITFUN with or without constrains), the fitted radiance is exactly equal to the radiance provided
;       There is probably something that I don't understand at this time concerning this problem. 
;       In any case, the first derivation of the temperature is probably not so bad (evaluation on going)  
;       
; EXAMPLES:
;
; HISTORY:
;	invtherm_H
;     SE, LESIA, Aug 2010:
;          Adapted from invSteins2
;	invtherm_H2
;     SE, LESIA, Sept 2014: adapted to 67P, merged spectra
;	invtherm_H2_mp & H3_mp
;     SE, LESIA, Sept 2014: adapted to 67P, merged spectra, use Mpfitfun (plots from v_invtherm4)
;      	Oct 2014: preserve input
;      	Feb 2016: Various fixes on the first step - output of temp (was always wrongly = input),
;				  fixed uncertainty (must always be >0), quiet option won't stop, etc
;				  added option/level kw
; Modification by B. Rousseau - 05/2018
;    -adaptation for VIRTIS-M, add options (rehabilitation of Powell function, cleaning the fitting functions, etc etc) 


; ******** Init *******
!EXCEPT=0     ; filters math error messages
if ~keyword_set(quiet) then  quiet= 0
if ~keyword_set(Lambert) then  Lambert = 0
If ~keyword_set(level) then  level= 1
If ~keyword_set(direct) then  direct=0 
If ~keyword_set(powel) then  powel=0
If ~keyword_set(delta) then  delta=0


; ******** sets output file names *******
;infile = 'ThermalHMP'
;outfile = getenv('Temp')+path_sep()+(strsplit(infile, '.' ,/ex))(0)+'_inv.idl'
;tempo = file_search(outfile)
;tt = 1
;While tempo do begin&$
; outfile =strcompress( (strsplit(infile, '.' ,/ex))(0) +'_inv' +string(tt) + '.idl', /re)&$
; tt=tt+1&$
; tempo = file_search(outfile)&$
;endwhile
;figname =strcompress('InverseRadHMP' +string(tt-1) + '.jpg', /re)

; Attributing the angles
If N_elements(m) EQ 0 then begin
  message, "No angles provided"
endif else begin
  mu= m
  mu0 = m0
  if ~quiet then print, "cosines:", mu0, mu
endelse

;Distance au Soleil et mise à l'échelle du flux solaire selon cette distance
if ~quiet then print, "Sun distance:", dista, ' au'
distS = dista
Rsoleil = Sflux
flux1= Rsoleil/(dista)^2 / !pi   ; Sun irradiance at target distance, correct scaling


; get rid of NaNs
Radiance = input_radiance ; input_radiance et radiance c'est la radiance du spectre en cours de traitement
ind = where(finite(Radiance) EQ 0, count) ; si jamais y a des NaN dans la radiance (? pas sur)
if count GT 0 then Radiance(ind) = 0.


; *** Merge orders as in plot_h, from merge
if VH eq 1 then begin ; because it means that we are dealing with H cube
  merge = merge_cube()
  if keyword_set(s_channel) then can = MERGE.S_TOOLS.S_CHANNEL
  if keyword_set(b_channel) then can = MERGE.B_TOOLS.MY_CHANNEL
  if ~keyword_set(b_channel) and ~keyword_set(s_channel) then begin
    print,'Channel keyword is missing, S_channel has been activated'
    can = MERGE.S_TOOLS.S_CHANNEL
  endif
  Radiance = Radiance(can)
  lam = lam(can)
  flux1 = flux1(can)
  Rsoleil =Rsoleil(can)
endif

;; smoothing & sampling - only for inversion/fit
;Radiance = (smooth(Radiance, 5))(1:*:5)
;lam = (smooth(lam, 5))(1:*:5)
;flux1 = (smooth(flux1, 5))(1:*:5)
;Rsoleil = (smooth(Rsoleil, 5))(1:*:5)	; used in commons!

; ******** simu rapide courts lam *******
; on définit la reflectance
Ntot = N_elements(lam)	; should be 1739 or 1645
refl = Radiance / (flux1) ; on divise la radiance (le spectre), par le flux solaire à la distance de l'acquisition, on a donc la réflectance
surfa = refl
lamS = lam
uncert = 0.05*radiance; incertitude donner comme 5%
deltaR = abs(radiance(1:Ntot-1)-radiance(0:Ntot-2))

; on définit une réflectance moyenne à partir du flux solaire à la bonne distance et en "scalant" sur la réflectance moyenne du spectre VIRTIS
; dans les longueurs d'onde entre 1.8 et 3.1µm 
; ça permet ensuite de dériver une première température à partir de cette réflectance estimée, basée sur une moyenne
; à vérif avec des température plus élevé, pas sur que ce soit good, mais a priori OK en vérifiant vite fait
if VH eq 1 then begin
  lamref1 = min(where(lam LE 2.5)) 	; lambda max où on peut considéré qu'il n'y a pas d'émission, OK
  If ~keyword_set(refmoy) then refmoy = mean(refl(lamref1:*))
endif
if VM eq 1 then begin
  lamref1 = max(where(lam LE 2.5))
  If ~keyword_set(refmoy) then refmoy = mean(refl(0:lamref1))
endif


uncert(where(uncert LE 0.)) = 0.05*refmoy	; là où l'incertitude pouvait être inf à 0 on la donne égale à 5% de la refl moyenne
 
 
if ~quiet then begin
; device, dec=0
; loadct, 12, /sil
; plot, lamS, refl*flux1, th=2, /xst, xtit= 'Wavelength ($\mu$m)', ytit = 'Radiance (SI units)', tit='Signal'
; oplot, lamS, flux1 * refmoy, col=100 ; estim réfléchi (en radiance)
 p = plot(lamS,refl*flux1, thick=2, xtit= 'Wavelength ($\mu$m)', ytit = 'Radiance (SI units)')
 p = plot(lamS, flux1*refmoy,color='dodger blue',/over, name='Flux sol. réfléchi approx')
 p.close
endif



; STEP 1
;**********************************************************************
;********* First fit for Temperature / from average reflectance *******
;**********************************************************************
; Dans cette première étape on évalue une première fois la température
; Ça permet ensuite d'évaluer au mieux le flux réfléchi du spectre en radiance pour l'enlever ensuite (STEP 2)
; et enfin de mieux ré-estimer la température (STEP 3)

; On définit selon l'option lambert, la fonction qui servira au fit de la réflectance
If Lambert eq 1 then fctname = 'inv1_Lambert' else fctname = 'inv1_LS' ; by default
if ~quiet then print, 'MPfitFun : '+fctname

; valeur initialisation
init=dblarr(2)
If keyword_set(temp) then init(0)=double(temp) else init(0) = 200.d ; estim T for 67P, soit elle est fournie, soiton prend 200K
init(1) = refmoy ; l'estimation du facteur de radiance pour démarrer est définie comme la moyenne de la réflectance entre 1.8 et 3.1 µm
 
if ~quiet then print, 'Initial guess of the temperature :',init(0)
if ~quiet then print, 'Initial guess of the reflectance (ref moy) : ',init(1)
ftol= 1.E-6 ; tolérance sur le fit


T0 = systime(1)
; T and mean refl only
if ~powel then begin
  retry:
  parms = MPFITFUN(fctname, LamS, Radiance, uncert, init, Niter=Niter1, YFIT= YFIT, BESTNORM= BESTNORM, PERROR= PERROR, STATUS= STATUS,ERRMSG=errmsg,quiet=1, MAXITER=500, Ftol=Ftol,DOF=DOF)
  if size(parms,/dim) eq 0 then begin
    print, 'Step 1 : ERROR'
    print, ERRMSG
    print, 'RETRY'
    goto, retry
  endif
  fun = Yfit
  Chi2= sqrt( total( (Radiance - Fun)^2 / uncert^2))
endif else begin
  xi = TRANSPOSE([[1.0, 0.0],[0.0, 1.0]])
  POWELL, init, xi, ftol, fmin, 'refinv1'  , itmax=50, iter=Niter1, /double     ; Lambert
  print, 'Parms via POWELL:', init  
  parms=init
endelse
if ~quiet then print, 'Fitted temperature 1 :', parms(0)


if ~quiet then begin
  print, 'Tp exec, 1 (s): ', systime(1) -T0
  print, 'Iter #:',niter1
  if ~powel then print, 'RMS radiance only: ', sqrt( total( (Radiance - Fun)^2 )/ Ntot)	
  if ~powel then print, 'Chi2 radiance only: ', chi2
  p = plot(lams,radiance,name='Radiance (data)')
  p = plot(lams,parms(1)*flux1, color='blue', /over,name='Flux sol. réfléchi estimé') ; plot du flux solaire à la distance de l'observation multiplié par la moyenne de la réflectance
  p = plot(lams,fun,/over,color='green',name='Fitted radiance')
  print, 'Fitted temperature:', parms(0)     
  print, 'Mean reflectance:', parms(1)     
;  DOF     = N_ELEMENTS(LamS) - N_ELEMENTS(PARMS) ; degrée de liberté
  PCERROR = PERROR * SQRT(BESTNORM / DOF) ;PERROR = 1-sigma pour chaque paramètre ; l'erreur calibré en fonction du problème, c'est à dir en fonction du numbre de degré de liberté. Voir doc mpfit
  print, 'Errors:', Pcerror     
endif

; On recalcule les grandeurs après le fit
Bt = corpsn(parms(0), lamS) ; W/sr/m2/µm - spectre de corps noir correspondant à la température dérivée
if lambert then emiss = 1. - (parms(1)) / mu0  else $ ; émissivité estimée d'après la loi de Kirchoff et d'après la réflectance moyenne corrigée du cosinus de l'incidence pour Lambert ou aussi de ceux des autres angles pour LS
  emiss = 1. - 2* parms(1) *(mu+mu0) /mu0 * (1. - mu*alog( (mu + 1.) / mu ))
; p=plot(lams,bt*emiss,/over,color='red')
; Repérage de la position où l'émission thermique n'est plus négligeable selon la première estimation de la température dérivée
if VH eq 1 then lamref = max(where(Bt*emiss/(parms(1)*flux1) GT 0.01))
if VM eq 1 then lamref = min(where(Bt*emiss/(parms(1)*flux1) GT 0.01))
if lamref EQ -1 then lamref = 1	; security in case of _very_ low T, on prendrait alors la première position

; To return
T1=parms(0)
rf1=parms(1)
if ~powel then chi2_1 = chi2 else chi2_1="POWEL"

If level EQ 1 then goto, suite; return T + fitted avg reflectance (refined from initial estimate)
  	; skip second level for the time being

; ************************************************************************************************************************************************
; ************************************************************************************************************************************************
; ************************************************************************************************************************************************


; STEP 2
;**********************************************************************
;************** Reflectance fit from Temperature estimate *************
;**********************************************************************
; ici on va estimer le facteur de radiance à partir du flux thermique que l'on calcule avec la température du fit précédent
; pour ça on retire le flux thermique à la radiance puis on divise par le flux solaire à la distance de 67P duquel on a retiré aussi le flux thermique, le tout corrigé des angles selon Lambert ou LS
; on part en fait de la formule de la radiance dans laquel on remplace l'émissivité par le facteur de radiance et à partir de la loi de Kirchoff

if ~direct then begin
  ; Méthode n° 1 : On utilise MPFITFUN aussi pour fitter la radiance et estimer la réflectance
  ; ça prend du temps et ça n'apporte pas de précision par rapport à un calcul direct (voir méthode n°2)
  T0 = systime(1) ; start time for computation
  if lambert then fctname = 'inv2_Lambert' else fctname = 'inv2_LS'
  init = replicate(parms(1), N_elements(Lams)) ; initialisation réflectance moyenne en fonction du fit précédent
  TempF = parms(0) ; initialisation température d'après fit précédent (common params6)
  print,"MPFITFUN used at the second step"

;*********************ON PREND QUE LE THERMIQUE POUR LE FIT ******************************
;  ; Lambda et radiance uniquement où le signal thermique n'est plus négligeable
  if VH eq 1 then begin
   lam_th = LamS(0:lamref)
   rad_th = Radiance(0:lamref)
   uncert_th = uncert(0:lamref)
   init_th = init(0:lamref)
  endif
  if VM eq 1 then begin
   lam_th = LamS(lamref:-1)
   rad_th = Radiance(lamref:-1)
   uncert_th = uncert(lamref:-1)
   init_th = init(lamref:-1)
  endif
;  
;  
;  ; commentaire Stéphane : "this fit is useless: always identical to direct computation within machine error"
  parms2 = MPFITFUN(fctname, lam_th, rad_th, uncert_th, init_th, Niter=Niter2, YFIT= YFIT2, BESTNORM= BESTNORM, PERROR= PERROR, STATUS= STATUS,ERRMSG=errmsg,/quiet, MAXITER=500, Ftol=Ftol, DOF= DOF)
;*********************ON PREND QUE LE THERMIQUE POUR LE FIT ******************************


;*********************AUTRE ESSAI EN PRENANT TOUT LE RANGE DE LAMBDA ******************************
; Lambda et radiance uniquement où le signal thermique n'est plus négligeable
;parinfo à tester
; commentaire Stéphane : "this fit is useless: always identical to direct computation within machine error"
;parms2 = MPFITFUN(fctname, lams, radiance, uncert, init, Niter=Niter2, YFIT= YFIT2, BESTNORM= BESTNORM, PERROR= PERROR, STATUS= STATUS,ERRMSG=errmsg,/quiet, MAXITER=500, Ftol=Ftol, DOF= DOF)
;*********************AUTRE ESSAI EN PRENANT TOUT LE RANGE DE LAMBDA ******************************

  fun = Yfit2
  Chi2= sqrt( total( (Radiance(0:lamref) - Fun)^2 / uncert(0:lamref)^2))  ; **** ça c'est louche ? ****
  if ~quiet then begin
    print, '  -- Level 2'
    print, '  Tp exec, 1 (s): ', systime(1) -T0
    print, 'Iter #:',niter2
    print, '  RMS radiance only: ', sqrt( total( (Radiance(0:lamref) - Fun)^2 )/ Ntot)
    print, '  Chi2 radiance only: ', chi2   ; = sqrt( total( (Radiance - Fun)^2 / uncert^2)), OK
    p = plot(lams,fun)
    p = plot(lams,Yfit2, col=100 ) ; overplots fit
;    DOF     = N_ELEMENTS(LamS(0:lamref)) - N_ELEMENTS(PARMS2)  ; OK, la fct retourne bien DOF = 0
;    PCERROR = PERROR * SQRT(BESTNORM / DOF)
;    print, 'Avg refl error:', mean(Pcerror)
  endif
endif else begin
  ; Méthode n°2 : Estimation directe de la réflectance selon Lambert ou LS
  ; on retire le flux thermique estimé à la radiance et on divise par le flux solaire à la distance de l'observation duquel on retire le flux thermique divisé par le cos(i)
  estrefl = (radiance-Bt)/(Flux1-(Bt/mu0)) ; estimation à la Lambert, simple
  ;estrefl_LS = ; compliqué, voir note, intégrale nécessaire
  Yfit2 =  Estrefl ; on  définit yfit2 à partir de l'estimation de la réflectance
endelse

if VH then begin
  parms2= Yfit2(0:lamref)
  init3= [parms(0),parms2, refl(lamref+1:*)] ; ok for H    ; add newly fit reflectance + reflected part
endif
if VM then begin
  parms2= Yfit2(lamref:-1)
  init3= [parms(0), refl(0:lamref-1),parms2]
endif
; return fitted T + fitted reflectance
; put back also filtered channels? 

;init3(0) = parms(0)	; pass T from level 1

If level EQ 2 then goto, suite	; skip third level 
; ************************************************************************************************************************************************
; ************************************************************************************************************************************************
; ************************************************************************************************************************************************


; STEP 3

; ************************************************************************
; ****************  T derivation again from reflectance fit   ************
; ************************************************************************
If Lambert and ~delta then fctname = 'inv3_Lambert'
If delta then begin
  If Lambert then fctname = 'invConst_Lambert' else fctname = 'invConst_LS'
Endif

parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0.D]}, N_elements(init3))
; température
parinfo[0].value = init3(0)
parinfo[0].fixed = 0  ; T is the variable
parinfo[0].limited = [1,1]
parinfo[0].limits = [150D,500D]  
; facteur de réflectance qui servira au calcul de l'émissivité 
parinfo[1:*].value = init3(1:*)
parinfo[1:*].fixed = 1	; radiance factor is fixed now

; T and fit refl 
parms3 = MPFITFUN(fctname, LamS, Radiance, uncert, init3, Niter=Niter3, YFIT= YFIT3, BESTNORM= BESTNORM, PERROR= PERROR, STATUS= STATUS,ERRMSG=errmsg, MAXITER=500, Ftol=Ftol, PARINFO=parinfo, /quiet,DOF=DOF)

;to return
T3 = parms3(0)
rf3 = parms3(1:-1)
Chi2_3= sqrt( total( (Radiance - YFIT3)^2 / uncert^2))

if ~quiet then print, 'Fitted temperature 3:', parms3(0)     








suite: 
if ~quiet then print, '  Tp exec, Tot (s): ', systime(1) -T0

if ~quiet then begin
  p=plot(lams, radiance,name ='Radiance') ; data ou flux1*refl (reflectance x flux solaire)
  p=plot(lams, flux1 * refmoy, color='blue',/over,name='Flux solaire réfléchi approx') ; ref moy entre 2 et 3µm *  flux solaire à la bonne distance
  p=plot(lams, init(1:*)*flux1, col='red', /over,name='Fitted radiance factor * Solar flux') ; facteur de radiance*flux solaire

  p=plot(lams,radiance-(init(1:*)*flux1),/over,color='blue',name='Radiance-(Rf * Solar flux)') ; test batiste
  p = plot(lams,bt,/over,color="red",name='Bt - Flux thermique Step 2') ; test batiste
  legend = legend(font_size = 20, /auto_text_color,SAMPLE_WIDTH=0.05)
  ; calcul inverse avec T estim�e, pour Lambert
  refEst = (Radiance - Bt) / (Flux1 - Bt / Mu0) ; calcul du facteur de radiance avec la température estimée
  ; c'est bien égal à init(1:*), cad le facteur de radiance en sorti d'étape 2
  p1=plot(lam, refl, thick = 2,name='Reflectance')
  p1=plot(lams, refEst, color = 'green',/over,name='Radiance factor Step 2') ; facteur de radiance recalculé
  p1=plot(lam, replicate(parms(1),N_elements(lams)), color='orange',/over,name='Constant radiance factor') ; la première estimation du facteur de radiance, celle qui est constante
  
  ; fit avec refl cst - permet de voir le CN
  reconstb =  refmoy +  corpsn(parms(0), lamS)   * (1.-refmoy/mu0) / ( Flux1 ) ; Bt reconstruite à partir de la première température
  p1=plot(lams, reconstb,col='red',/over,name='Rebuild thermal flux T from STEP 1')
  ; estim correcte
  
  emiss = 1. - refEst / mu0
  p1=plot(lams, emiss, color = 'violet',/over,name='Emissivity STEP 2')
  ; Reconstruction OK avec mod�le Lambert 
  
  reconst =  (refEst * Flux1 +  Bt * emiss) / ( Flux1 ) ; réflectance
  p1=plot(lams, reconst, color='blue',/over,name='Reflectance STEP 2') ; même chose que refl

  ; avec température initiale
  Talt = parms(0)
  reconst2 =  refmoy +  corpsn(Talt, lamS)   *  (1.-refmoy/mu0) / ( Flux1 )
  p1=plot(lams, reconst2, color = 'blue' , linest=2,/over,name='Rebuild reflectance') ; reflectance reconstruite
  refEst2 = (Radiance - corpsn(Talt, lamS)) / (Flux1 - corpsn(Talt, lamS)  / Mu0)
  p1=plot(lams, refest2, color = 'violet' , linest=2,/over, name = 'Radiance factor') ; facteur de radiance
  p1=plot(lams, (1.-refEst2/mu0), color = 'green', linest=2,/over, name = 'Emissivity from 1st T') ; émissivité
  legend = legend(font_size = 20, /auto_text_color,SAMPLE_WIDTH=0.05)
endif

END 





; ******************************************************************************************************************************
; ******************************************************************************************************************************
; ******************************************************************************************************************************





; ****************** FONCTION DE LAMBERT et LOMMEL-SEELIGER POUR LES FIT ******************
; FIRST ESTIMATION OF THE REFLECTANCE   
;STEP 1
FUNCTION inv1_Lambert, lamS, x ; Estimation à la Lambert
  common param5, fluxS, distS, mu, mu0
  
  T= x(0)
  surfa = replicate(x(1),N_elements(lams))     ; radiance factor
  Bt = corpsn(T, lamS)     ; W / sr / m2 / �m
  emiss = 1. - surfa / mu0     ; pas fct des angles (cst avec mu0)
  fun =  surfa * FluxS/!pi/distS^2 +  Bt * emiss
  return, fun
end

FUNCTION inv1_LS, lamS, x ; estimation avec Lommel-Seeliger
  common param5, fluxS, distS, mu, mu0
  
  T= x(0) ; température initial, guess
  surfa = replicate(x(1),N_elements(lams)) ; radiance factor => constant basé sur la réflectance moyenne
  Bt = corpsn(T, lamS) ; W / sr / m2 /  spectre de corps noir en radiance à partir de l'estimation de la température
  emiss = 1. - 2* surfa *(mu+mu0) /mu0 * (1. - mu*alog( (mu + 1.) / mu )) ; calcul d'une émissitivé constante à partir d'un radiance factore constant et des angles (modèle de Lommel Seeliger ici) 
  fun =  surfa * FluxS/!pi/distS^2 +  Bt * emiss ; simulation de la radiance basée sur la réflectance moyenne calculé à partir de la moyenne de la radiance mesuré à courte lambda (hyp : pas de flux thermique ici), de l'émissivité basée sur cette réflectance et enfin sur le spectre du corps noir avec la températur fournie en hyp
  return, fun
end

; fonction utilisé par POWELL
FUNCTION refinv1, x
  common param2, flux1, lams
  common res, fun
  common param5, RSoleil, distS, mu, mu0
  common param7, Radiance

  ; estimation a la Lambert
  Ntot = N_elements(lamS)
  T= x(0)
  surfa = x(1:*) ; refl directionnelle
  Bt = corpsn(T, lamS)     ; W / sr / m2 / �m
  emiss = 1. - surfa / mu0  ; pas fct des angles (cst avec mu0)
  fun =  surfa * Flux1 +  Bt * emiss
  
  return, sqrt( total( (Radiance - Fun)^2 )/ Ntot); on retourne le chi², ce qu'on essaie de minimiser
end



;STEP 2
FUNCTION inv2_Lambert, lamS, x  ; Lambert
  common param5, fluxS, distS, mu, mu0
  common param6, TempF
  common param7, Radiance

  Bt = corpsn(TempF, lamS)     ; W / sr / m2 / �m
  emiss = 1. - x / mu0     ; pas fct des angles (cst avec mu0), on calcul l'émissivité à partir du facteur de réflectance qui va bouger selon les itérations
  fun =  x * FluxS/!pi/dists^2 +  Bt * emiss ; on calcul la radiance puis le programme estime la différence avec les données etc etc jusqu'au meilleir fit
;  p=plot(lams,radiance)
;  p=plot(lams,fun,/over,color='red')
;  p=plot(lams,1-emiss,/over,color='blue')
;  p.close
  return, fun
end


FUNCTION inv2_LS, lamS, x  ; L-S
  common param5, fluxS, distS, mu, mu0
  common param6, TempF

  ; estimation style Lommel-Seeliger
  ; Ntot =N_elements(lamS)
  ;T= x(0)
  ;surfa = replicate(x(1),N_elements(lams))     ; radiance factor
  Bt = corpsn(TempF, lamS)     ; W / sr / m2 / �m
  emiss = 1. - 2* x *(mu+mu0) /mu0 * (1. - mu*alog( (mu + 1.) / mu ) )
  fun =  x * FluxS/!pi/distS^2 +  Bt * emiss
  ; oplot, lamS, fun
  ;return, sqrt( total( (Radiance - Fun)^2 )/ Ntot)
  return, fun
end






FUNCTION inv3_Lambert, lamS, x ; Lambert
  common param5, fluxS, distS, mu, mu0
  common param6, TempF
  common param7, Radiance

  ; estimation a la Lambert
  ;surfa = x     ; radiance factor
  Bt = corpsn(x(0), lamS)     ; W / sr / m2 / �m
  emiss = 1. - x(1:*) / mu0     ; pas fct des angles (cst avec mu0)
  fun =  x(1:*) * FluxS/!pi/distS^2 +  Bt * emiss
;  p=plot(lams,radiance,thick=2)
;  p=plot(lams,fun,/over,color="red")
;  print, sqrt( total( (Radiance - Fun)^2 / uncert^2)) ; attention aux limites / uncert
;  p.close
  return, fun
end









FUNCTION invConst_LS, lamS, x    ; L-S avec contrainte
  ; x = [T, refl]     variables ind�pendantes
  common param1, uncert, deltaR, Lcrack
  common param5, fluxS, distS, mu, mu0
  common param7, Radiance
  common param2, flux1
  common res1, fun

  Ntot =N_elements(lamS)
  T= x(0)
  surfa = x(1:Ntot)     ; radiance factor
  Bt = corpsn(T, lamS)     ; W / sr / m2 / �m
  emiss = 1. - 2* surfa *(mu+mu0) /mu0 * (1. - mu*alog( (mu + 1.) / mu ) ) ; LS

  fun =  surfa * FluxS/!pi/distS^2 +  Bt * emiss
  ;delta = x(Ntot+1:*)
;  delta1 = abs(x(2:Ntot)-x(1:Ntot-1))
  delta1 = abs(surfa(1:Ntot-1)-surfa(0:Ntot-2)) ; définit comme ligne au dessus mais à partir de surfa, donc même syntaxe que deltaR
  fun = [fun,delta1(Lcrack:*)]
  ;oplot, lamS, fun
  ;oplot, lamS, x(2:Ntot)-x(1:Ntot-1), col=100

  ;plot, lamS, radiance / (FluxS/!pi/distS^2), /xst
  ;oplot, lamS, surfa, col=30     ; refl surf ajust�e
  ;oplot, lamS, deltaR, col=100
  ;oplot, lamS, uncert(Ntot:*), col=80
  ;oplot, lamS, delta1, col=200     ; delta refl surf ajust�e (doit �tre bleu � cyan)

  ;print, sqrt( total( (Radiance - Fun(0:Ntot-1))^2 / uncert(0:Ntot-1)^2))
  ;print, sqrt( total( (deltaR - delta1)^2 / uncert(Ntot:*)^2))
  ;return, sqrt( total( (Radiance - Fun)^2 )/ Ntot)
  ; return, sqrt( total( ([Radiance,deltaR(Lcrack:*)] - Fun)^2 / uncert^2))
  return, sqrt( total( ([Radiance,deltaR] - Fun)^2 / uncert^2))     ; TBC, mais �a devrait aller�
end


FUNCTION invConst_Lambert, lamS, x    ; Lambert avec contrainte
  ; x = [T, refl]     variables indépendantes
  common param1, uncert, deltaR, Lcrack
  common param5, fluxS, distS, mu, mu0
  common param7, Radiance
  common param2, flux1
  common res1, fun

  Ntot =N_elements(lamS)
  T= x(0)
  surfa = x(1:Ntot)     ; radiance factor
  Bt = corpsn(T, lamS)     ; W / sr / m2 / �m
  emiss = 1. - surfa / mu0     ; pas fct des angles (cst avec mu0) - Lambert
  
  fun =  surfa * FluxS/!pi/distS^2 +  Bt * emiss
  
;  p=plot(lams,bt,/over,color='red')
  
  ;delta = x(Ntot+1:*) ; c'est quoi ça ? un premier test ? 
;  delta1 = abs(x(2:Ntot)-x(1:Ntot-1))
  delta1 = abs(surfa(1:Ntot-1)-surfa(0:Ntot-2)) ; définit comme la ligne au dessus mais à partir de surfa, donc même syntaxe que deltaR
;  fun = [fun,delta1(Lcrack:*)] ; qu'est-ce que Lcrack ? 
  fun = [fun,delta1]
;  print,sqrt( total( ([Radiance,deltaR] - Fun)^2 / uncert^2))
  return, sqrt( total( ([Radiance,deltaR] - Fun)^2 / uncert^2))  ; qu'est-ce que deltaR, la différence entre la premier fit et les données ?    ; TBC, mais �a devrait aller�
end