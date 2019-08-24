pro weird_backups

 
  
  DIRDATA='/Users/fandrieu/data/'
  STP=67
  ;num='H1_00396778317'
  dyn_clip=1
  display_clip=1
  save_p=0
  ratio=1
  list_H_cube=['H1_00396778317', 'H1_00396818786', 'H1_00397037997', $
    'H1_00397077249', 'H1_00397243610', 'H1_00397294310', 'H1_00396782306', $
    'H1_00396825340', 'H1_00397041036', 'H1_00397083369', 'H1_00397258490', $
    'H1_00397300250', 'H1_00396788782', 'H1_00396837576', 'H1_00397047178', $
    'H1_00397095314', 'H1_00397264070', 'H1_00397312610', 'H1_00396800333', $
    'H1_00396841325', 'H1_00397059123', 'H1_00397138801', 'H1_00397276430' ,$
    'H1_00397316997', 'H1_00396804897', 'H1_00396855961', 'H1_00397062874', $
    'H1_00397147950', 'H1_00397279790', 'H1_00397330046']
  
  num=list_H_cube[0]
  dark_number=1
  
  
  restore, '/Users/fandrieu/Documents/Programmes/IDL_sav_files/CT_custom.sav'

  cube_name=strcompress(DIRDATA+'/STP'+string(stp, format='(I03)')+$
    '/'+num+'.QUB', /remove_all)
  cube_name_cal=strcompress(DIRDATA+'/STP'+string(stp, format='(I03)')+$
    '/'+num+'.CAL', /remove_all)
  cube_name_dark=strcompress(DIRDATA+'/STP'+string(stp, format='(I03)')+$
    '/'+num+'.DRK', /remove_all)
    
  cube_raw = virtispds(cube_name, /silent)
  cube_cal = virtispds(cube_name_cal, /silent)
  cube_dark_cal = virtispds(cube_name_dark, /silent)
 
 
  
  ;ploplo = plot_virtish(spectre_str_a.data[*,6])
  ploplo = plot_virtish(cube_cal.qube[*,5*dark_number+1])
  ploplo = plot_virtish(cube_dark_cal.qube[*,5*dark_number+1])
  ;tvscl,cube_raw.qube[*,*,5]
  ;tvscl,congrid(cube_raw.qube[0:430:2,140:255,5]-median(cube_raw.qube[0:430:2,140:255,5]),700, 300)<200>0
  ;tvscl,congrid(cube_raw.qube[0:430:2,0:140,5]-median(cube_raw.qube[0:430:2,0:140,5]),700, 200)>0<1000
  

  wawa=cube_raw.qube[0:430:2,0:140,5*dark_number]
  sowa=sort(wawa)
  n_el=n_elements(wawa)
  n_clip=fix(n_el*dyn_clip*0.01)

  ;mimi=median(wawa[sowa[0:n_clip]])
  ;mama=median(wawa[sowa[n_el-n_clip:n_el-1]])
  mimi=wawa[sowa[n_clip]]
  mama=wawa[sowa[n_el-n_clip-1]]
  if (display_clip ne 0) then begin
    wawa[where(wawa le mimi)]=mimi
    wawa[where(wawa ge mama)]=mama
  endif
  fullrange=mama-mimi
  mimi=mimi-fullrange*0.01
  mama=mama+fullrange*0.01
  c1 = CONTOUR(wawa, /FILL, RGB_TABLE=ct, aspect_ratio=0.5*ratio,n_levels=25, $
    TITLE=num+' DARK ', POSITION=[0.05,0.18,0.95,0.9], min_value=mimi, max_value=mama)
  cb1 = COLORBAR(TITLE='DARK',TARGET=c1, ORIENTATION=0, $
    POSITION=[0.3,0.08,0.7,0.12], taper=1, TICKFORMAT='(F7.0)')

  if (save_p eq 1) then begin
    plotname=strcompress(dirplot+'/T1_'+num+'_CALIBRATED_'+'at_'+string(cube_cal.table[0,posw])+'.png', /remove_all)
    c1.Save, plotname, BORDER=10, RESOLUTION=300, /TRANSPARENT
  endif

  wawa=cube_raw.qube[0:430:2,140:255,5*dark_number]
  sowa=sort(wawa)
  n_el=n_elements(wawa)
  n_clip=fix(n_el*dyn_clip*0.01)

  ;mimi=median(wawa[sowa[0:n_clip]])
  ;mama=median(wawa[sowa[n_el-n_clip:n_el-1]])
  mimi=wawa[sowa[n_clip]]
  mama=wawa[sowa[n_el-n_clip-1]]
  if (display_clip ne 0) then begin
    wawa[where(wawa le mimi)]=mimi
    wawa[where(wawa ge mama)]=mama
  endif
  fullrange=mama-mimi
  mimi=mimi-fullrange*0.01
  mama=mama+fullrange*0.01
  c1 = CONTOUR(wawa, /FILL, RGB_TABLE=ct, aspect_ratio=0.5*ratio,n_levels=25, $
    TITLE=num+' DARK ', POSITION=[0.05,0.18,0.95,0.9], min_value=mimi, max_value=mama)
  cb1 = COLORBAR(TITLE='DARK',TARGET=c1, ORIENTATION=0, $
    POSITION=[0.3,0.08,0.7,0.12], taper=1, TICKFORMAT='(F7.0)')


;
;  wawa=cube_raw.qube[0:430:2,*,5*dark_number]
;  sowa=sort(wawa)
;  n_el=n_elements(wawa)
;  n_clip=fix(n_el*dyn_clip*0.01)
;
;  ;mimi=median(wawa[sowa[0:n_clip]])
;  ;mama=median(wawa[sowa[n_el-n_clip:n_el-1]])
;  mimi=wawa[sowa[n_clip]]
;  mama=wawa[sowa[n_el-n_clip-1]]
;  if (display_clip ne 0) then begin
;    wawa[where(wawa le mimi)]=mimi
;    wawa[where(wawa ge mama)]=mama
;  endif
;  fullrange=mama-mimi
;  mimi=mimi-fullrange*0.01
;  mama=mama+fullrange*0.01
;  c1 = CONTOUR(wawa, /FILL, RGB_TABLE=ct, aspect_ratio=0.5*ratio,n_levels=25, $
;    TITLE=num+' DARK ', POSITION=[0.05,0.18,0.95,0.9], min_value=mimi, max_value=mama)
;  cb1 = COLORBAR(TITLE='DARK',TARGET=c1, ORIENTATION=0, $
;    POSITION=[0.3,0.08,0.7,0.12], taper=1, TICKFORMAT='(F7.0)')
;
;  wawa=cube_raw.qube[*,*,5*dark_number]
;  sowa=sort(wawa)
;  n_el=n_elements(wawa)
;  n_clip=fix(n_el*dyn_clip*0.01)
;
;  ;mimi=median(wawa[sowa[0:n_clip]])
;  ;mama=median(wawa[sowa[n_el-n_clip:n_el-1]])
;  mimi=wawa[sowa[n_clip]]
;  mama=wawa[sowa[n_el-n_clip-1]]
;  if (display_clip ne 0) then begin
;    wawa[where(wawa le mimi)]=mimi
;    wawa[where(wawa ge mama)]=mama
;  endif
;  fullrange=mama-mimi
;  mimi=mimi-fullrange*0.01
;  mama=mama+fullrange*0.01
;  c1 = CONTOUR(wawa, /FILL, RGB_TABLE=ct, aspect_ratio=ratio, n_levels=25,$
;    TITLE=num+' DARK ', POSITION=[0.05,0.18,0.95,0.9], min_value=mimi, max_value=mama)
;  cb1 = COLORBAR(TITLE='DARK',TARGET=c1, ORIENTATION=0, $
;    POSITION=[0.3,0.08,0.7,0.12], taper=1, TICKFORMAT='(F7.0)')


end