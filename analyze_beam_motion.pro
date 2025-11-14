; This procedure reads and HDF5 file containing a time-series of images of the X-ray beam on a scintillator
; It does the following:
; 1) Reads the specified HDF5 file 
; 2) Selects an ROI to analyze, defined by xrange and yrange
; 3) Determines the time per image from the timestamps of the first and last image in the file.
; 4) For each image:
;   - Averages all of the columns in the image to get an average vertical profile
;   - Computes the center of gravity (COG) of the vertical profile
; 5) Subtracts the mean COG, so the COG varies around 0
; 6) Plots the center of gravity as a function of time, annotated with the peak-to-peak and RMS variation
; 7) Computes the FFT of the time-series of the COG 
; 8) Plots the power-spectrum of the FFT
; 9) Repeats steps 4-8 for the horizontal direction
pro analyze_beam_motion, file, title, xrange, yrange, pixel_size

  data_full = h5_getdata(file, '/exchange/data')
  ; Extract ROI near the beam
  image_display, data_full[*,*,0]
  data = data_full[xrange[0]:xrange[1], yrange[0]:yrange[1], *]
  s = size(data, /dimensions)
  nx = s[0]
  ny = s[1]
  nimages = s[2]
  ;unique_id = h5_getdata(file, '/entry/instrument/NDAttributes/NDArrayUniqueId')
  ;time = h5_getdata(file, '/entry/instrument/NDAttributes/NDArrayTimeStamp')
  unique_id = h5_getdata(file, '/defaults/NDArrayUniqueId')
  time = h5_getdata(file, '/defaults/NDArrayTimeStamp')
  time_step = (time[nimages-1]-time[0]) / (nimages-1)
  time_axis = findgen(nimages) * time_step
  sampling_freq = 1./time_step
  nyquist_freq = sampling_freq/2.
  nfreq = nimages/2
  freq_axis = findgen(nfreq)/(nfreq) * nyquist_freq
  

  loadct, 5
  ; fps and bps seem to have no effect on these low-res movies
  ;make_movie, data, mp4_file=title + '_movie.mp4', /color, fps=sample_freq, bps=1e6
  make_movie, data, mp4_file=title + '_movie.mp4', /color
  
  ; Compute vertical beam motion
  ave = mean(data, dimension=1, /double)
  cog = dblarr(nimages)
  x = findgen(nimages)
  for i=0, nimages-1 do begin
    col = ave[*,i]
    cog[i] = total(col*x) / total(col)
  endfor
  ; Subtract mean and convert from pixels to microns
  cog = (cog-mean(cog)) * pixel_size
  peak_to_peak = max(cog)-min(cog)
  rms = stddev(cog)
  fmt = '(F6.2)'
  p = plot(time_axis, cog, xtitle='Time (s)', ytitle='Vertical beam position microns (center-of-gravity)', $
           title='Vertical beam position in ' + title + '(' + file + ')')
  t = text(.2, .8, 'Peak to peak= ' + string(peak_to_peak,format=fmt) + ' microns,   RMS= ' + string(rms, format=fmt) + ' microns')
  p.save, title + '_vertical_beam_motion_COG.png'
  
  power_spectrum = abs(fft(cog))
  ; Only use first half
  power_spectrum = power_spectrum[0:nfreq-1]
  ; Remove 0 frequency
  power_spectrum[0] = 0
  p = plot(freq_axis, power_spectrum, xtitle='Frequency (Hz)', ytitle='Vertical power spectrum', $
    title='Vertical power spectrum in ' + title + '(' + file + ')')
  p.save, title + '_vertical beam_motion_PS.png'
  
  ; Compute horizontal beam motion
  ave = reform(mean(data, dimension=2, /double))
  cog = dblarr(nimages)
  x = findgen(nimages)
  for i=0, nimages-1 do begin
    row = ave[*,i]
    cog[i] = total(row*x) / total(row)
  endfor
  ; Subtract mean and convert from pixels to microns
  cog = (cog-mean(cog)) * pixel_size
  peak_to_peak = max(cog)-min(cog)
  rms = stddev(cog)
  fmt = '(F6.2)'
  p = plot(time_axis, cog, xtitle='Time (s)', ytitle='Horizontal beam position microns (center-of-gravity)', $
    title='Horizontal beam position in ' + title + '(' + file + ')')
  t = text(.2, .8, 'Peak to peak= ' + string(peak_to_peak,format=fmt) + ' microns,   RMS= ' + string(rms, format=fmt) + ' microns')
  p.save, title + '_horizontal beam_motion_COG.png'
   
  power_spectrum = abs(fft(cog))
  ; Only use first half
  power_spectrum = power_spectrum[0:nfreq-1]
  p = plot(freq_axis, power_spectrum, xtitle='Frequency (Hz)', ytitle='Horizontal power spectrum', $
    title='Horizontal beam motion in ' + title + '(' + file + ')')
  p.save, title + '_horizontal_beam_motion_PS.png'
  
end
