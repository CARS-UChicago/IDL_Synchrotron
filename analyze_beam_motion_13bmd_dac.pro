; Program to analyse beam motion at 13-BM-D DAC

;title = '13BMD_DAC_Feb18, 10Hz'
;file = 'stream_fast_mono8_10Hz_001.h5'
;title = '13BMD_DAC_Feb18, 80Hz'
;file = 'stream_fast_mono8_80fs_002.h5'
;title = '13BMD_DAC_Feb18, 150Hz'
;file = 'stream_fast_mono8_150fs_003.h5'
;title = '13BMD_DAC_Mar4, Chiller On, 150Hz'
;file = 'stream_fast_mono8_150fs_chiller-ON_005.h5'
title = '13BMD_DAC_Mar4, Chiller Off, 150Hz'
file = 'stream_fast_mono8_150fs_chiller-off_006.h5'
;xrange = [860, 1060]
;xrange = [980, 1180]
xrange = [1300, 1500]
yrange = [0,299]
pixel_size = 0.2
analyze_beam_motion, file, title, xrange, yrange, pixel_size
end
