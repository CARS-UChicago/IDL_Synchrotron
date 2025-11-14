; Program to analyse vertical motion at 13-ID-C after APS-U
; 
file = 'IDC_mono-39keV_002.h5'
data_full = h5_getdata(file, '/entry/data/data')
; Extract ROI near the beam
data = data_full[750:1800, 0:1199, *]
s = size(data, /dimensions)
nx = s[0]
ny = s[1]
nimages = s[2]
pixel_size = 1.09
unique_id = h5_getdata(file, '/entry/instrument/NDAttributes/NDArrayUniqueId')
time = h5_getdata(file, '/entry/instrument/NDAttributes/NDArrayTimeStamp')
time_step = (time[nimages-1]-time[0]) / (nimages-1)
time_axis = findgen(nimages) * time_step
sampling_freq = 1./time_step
nyquist_freq = sampling_freq/2.
nfreq = nimages/2
freq_axis = findgen(nfreq)/(nfreq) * nyquist_freq

p = image(data[*,*,0], min=0, max=100)
p.save, '13IDC_beam_motion_Image ('+file+').png'

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
p = plot(time_axis, cog, xtitle='Time (s)', ytitle='Vertical beam position microns (center-of-gravity)', $
         title='Vertical beam position in 13-ID-C ('+file+')')
p.save, '13IDC_vertical_beam_motion_COG_'+file+'.png'

power_spectrum = abs(fft(cog))
; Only use first half
power_spectrum = power_spectrum[0:nfreq-1]
; Remove 0 frequency
power_spectrum[0] = 0
p = plot(freq_axis, power_spectrum, xtitle='Frequency (Hz)', ytitle='Power spectrum', $
  title='Vertical power spectrum in 13-ID-C ('+file+')')
p.save, '13IDC_vertical beam_motion_'+file+'.png'

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
p = plot(time_axis, cog, xtitle='Time (s)', ytitle='Horizontal beam position microns (center-of-gravity)', $
  title='Horizontal beam position in 13-ID-C ('+file+')')
p.save, '13IDC_horizontal_beam_motion_COG_'+file+'.png'

power_spectrum = abs(fft(cog))
; Only use first half
power_spectrum = power_spectrum[0:nfreq-1]
p = plot(freq_axis, power_spectrum, xtitle='Frequency (Hz)', ytitle='Power spectrum', $
  title='Horizontal power spectrum in 13-ID-C ('+file+')')
p.save, '13IDC_horizontal_beam_motion_PS_'+file+'.png'

end
