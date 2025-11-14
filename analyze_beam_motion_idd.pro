; Program to analyse beam motion at 13-ID-D

;title = '13IDD_November_22_cryocooler_off'
;title = '13IDD_November_22_high-frame_rate'
title = '13IDD_Feb11_Mono16
;file = 'statistics_cryocooler_001.h5'
;file = 'statistics_high-frame_002.h5'
file = 'idd_pg2_11Feb25_001_Mono16.h5'
xrange = [925, 1025]
;yrange = [550, 650]
yrange = [0, 99]
pixel_size = 0.349
analyze_beam_motion, file, title, xrange, yrange, pixel_size
end
