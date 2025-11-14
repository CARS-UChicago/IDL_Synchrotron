; Program to analyse beam motion at 13-BM-D April 9, 2025

;title = '13BMD_White_Apr9, 164Hz'
;file = 'white_164fps_0p5mmCu.h5'
;title = '13BMD_Mono_Apr9, 164Hz'
;file = 'mono_164fps_28keV.h5'
title = '13BMD_Mono_Apr9, Water Off, 164Hz'
file = 'white_164fps_28keV_water_off.h5'
xrange = [0, 1919]
yrange = [0,1199]
pixel_size = 5.74
analyze_beam_motion, file, title, xrange, yrange, pixel_size
end
