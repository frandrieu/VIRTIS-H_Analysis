pro stat_virtis_h


  ;distant_server=0
  dirname='/Users/fandrieu/Documents/Programmes/IDL_sav_files'
  ;if  distant_server eq 1 then dirname=strcompress('fandrieu@lesia08.obspm.fr:'+dirname, /remove_all)
  order=8
  firstmtp=4
  lastmtp=35
  nbwl=3456
  nbmtp=lastmtp-firstmtp
  stat_mom_s=fltarr(nbmtp, nbwl,4)
  stat_mom_c=fltarr(nbmtp, nbwl,4)
  for mtp_number=firstmtp,lastmtp-1 do begin
    indice=mtp_number-firstmtp
    MTP_name=string( mtp_number, format='(I2.2)')
    savname=strcompress(dirname+'/save_moments-MTP0'+MTP_name+'.sav', /remove_all)
    restore, savname
    if n_elements( mo_s) gt 1 then stat_mom_s[indice, *,*]= mo_s else  stat_mom_s[indice, *,*]=!Values.F_NAN
    if n_elements( mo_c) gt 1 then stat_mom_c[indice, *,*]= mo_c else  stat_mom_c[indice, *,*]=!Values.F_NAN

  endfor
  window,0
  device, decomposed=0
  loadct,39
  erase,255
  !P.background=0
  plot= plot_virtish( stat_mom_c[0,*,0]/stat_mom_c[0,*,1])
  loadct,39
  if (order gt 7) then begin
    for mtp_number=firstmtp+1,lastmtp-1 do begin
      indice=nbmtp-mtp_number+firstmtp
      plot= plot_virtish( stat_mom_c[indice,*,0]/stat_mom_c[indice,*,1], /over,  color=fix(255.*float(mtp_number-firstmtp)/float(nbmtp))) 
    endfor
  endif else begin
    plot= plot_virtish( stat_mom_c[1,*,0]/stat_mom_c[1,*,1], order=order)
    for mtp_number=firstmtp+1,lastmtp-1 do begin
      indice=nbmtp-mtp_number+firstmtp
      plot= plot_virtish( stat_mom_c[indice,*,0]/stat_mom_c[indice,*,1],  order=order, /over,  color=fix(255.*float(mtp_number-firstmtp)/float(nbmtp)))
    endfor
  endelse
  breakbreak=pointpoint
end
