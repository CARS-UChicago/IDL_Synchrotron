; Program to analyse vertical motion at 13-BM-D after APS-U

infile = '8_bit_163fps.h5'
out_prefix = '13BMD_163fps'
data = h5_getdata(infile, '/exchange/data')
s = size(data, /dimensions)
nx = s[0]
ny = s[1]
nimages = s[2]
pixel_size = 5.74
unique_id = h5_getdata(infile, '/defaults/NDArrayUniqueId')
time = h5_getdata(infile, '/defaults/NDArrayTimeStamp')
time_step = (time[nimages-1]-time[0]) / (nimages-1)
time_axis = findgen(nimages) * time_step
sampling_freq = 1./time_step
nyquist_freq = sampling_freq/2.
nfreq = nimages/2
freq_axis = findgen(nfreq)/(nfreq) * nyquist_freq

p = image(data[*,*,0])
p.save, '13BMD_beam_motion_Image.png'
 
ave = mean(data, dimension=1, /double)
x_mm = findgen(ny)*pixel_size/1000.
p = plot(x_mm, ave[*,0], xrange=[0,7], xtitle='Vertical position (mm)', ytitle='Beam intensity', $
  title="Average vertical beam profiles (every 1000'th image)")
p1 = plot(x_mm, ave[*, 1000], color='blue', over=p)
p1 = plot(x_mm, ave[*, 2000], color='red', over=p)
p1 = plot(x_mm, ave[*, 3000], color='green', over=p)
p1 = plot(x_mm, ave[*, 4000], color='purple', over=p)
p1 = plot(x_mm, ave[*, 5000], color='magenta', over=p)
p1 = plot(x_mm, ave[*, 6000], color='cyan', over=p)
p1 = plot(x_mm, ave[*, 7000], color='pink', over=p)
p1 = plot(x_mm, ave[*, 8000], color='orange', over=p)
p1 = plot(x_mm, ave[*, 9000], color='lavender', over=p)
p1 = plot(x_mm, ave[*, 9999], color='yellow', over=p)
p.save, out_prefix+'_vertical profile.png'

cog = dblarr(nimages)
x = findgen(nimages)
for i=0, nimages-1 do begin
  col = ave[*,i]
  cog[i] = total(col*x) / total(col)
endfor
; Subtract mean and convert from pixels to microns
cog = (cog-mean(cog)) * pixel_size
print, 'Peak-to-peak variation in COG=', max(cog)-min(cog)
print, 'RMS variation in COG=', stddev(cog)
p = plot(time_axis, cog, xrange=[0,62], xtitle='Time (s)', ytitle='Beam position (center-of-gravity, microns)', $
         title='White beam position in 13-BM-D @58m')
p.save, out_prefix+'_beam_motion_COG.png'

p = plot(time_axis[0:162], cog[0:162], xtitle='Time (s)', ytitle='Beam position (center-of-gravity, microns)', $
  title='White beam position in 13-BM-D @58m', symbol='Circle', sym_size=.5)
p.save, out_prefix+'_beam_motion_COG_1sec.png'

power_spectrum = abs(fft(cog))
; Only use first half
power_spectrum = power_spectrum[0:nfreq-1]
p = plot(freq_axis, power_spectrum, xrange=[0,85], xtitle='Frequency (Hz)', ytitle='Power spectrum', $
  title='White beam power spectrum @58m')
p.save, out_prefix+'_beam_motion_PS.png'

big = where(power_spectrum gt 1)
print, 'Power spectrum frequencies with amplitude greater than 1'
for i=0,n_elements(big)-1 do print, freq_axis[big[i]], power_spectrum[big[i]]

end
