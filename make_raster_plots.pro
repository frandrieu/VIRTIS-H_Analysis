pro make_raster_plots

  restore, 'Documents/Programmes/IDL_sav_files/raster_list_dominique.sav'
  wavelength=3.1
  n_begin=1
  n_end=n_elements(raster_list)-1
  n_end=1
  
  for rastn=n_begin, n_end do plot_XY_raster_VirtisH, raster_list[rastn], $
    wavelength, order=2, dyn_clip=0.5, /display_clip, /save_p, $
     dirname='/Users/fandrieu/Documents/Programmes/IDL_sav_files/Plots/RASTERS/Wave4a'
  for rastn=n_begin, n_end do plot_XY_raster_VirtisH, raster_list[rastn], $
    wavelength, order=3, dyn_clip=0.5, /display_clip, /save_p, $
    dirname='/Users/fandrieu/Documents/Programmes/IDL_sav_files/Plots/RASTERS/Wave4b'

end