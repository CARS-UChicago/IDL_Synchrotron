pro shielding_calcs, file

; Read and parse the list of parameters which the HTML form generated
openr, 1, file+".inp"
line = ""
readf, 1, line
close, 1
n = strparse(line, "=&", words)
n = n_elements(words)
words = reform(words, 2, n/2)
type = words(1,0)
pb = float(words(1,1))
xdist = float(words(1,2))
zend = float(words(1,3))
zdist = float(words(1,4))
xend = float(words(1,5))

; Set up color tables
black = 0
white = 1
red = 2
green = 3
blue = 4
yellow = 5
red_lut = [0, 255, 255, 0, 0, 255]
green_lut = [0, 255, 0, 255, 0, 255]
blue_lut = [0, 255, 0, 0, 255, 0]

; We go through the following twice, plotting once to GIF and once to PS
for i=0, 1 do begin
  if (i eq 0) then begin
     set_plot, 'ps'
     device, xsize=10, ysize=7, /inch, /land, /color, bits=4, file=file+".ps"
     tvlct, red_lut, green_lut, blue_lut
     !p.font=0
     !p.thick=3
     !p.charsize=1.5
     chsize=1.5
  endif else begin
     set_plot, 'z'
     !p.font=-1
     !p.charsize=1.5
     !p.thick=1
     chsize=1.5
  endelse

  if (type eq 'back') then begin
     X = 0.1 + xend*findgen(100)/99.
     brem_back, zdist, X, pb, Thick, Thin_target
  
     ymax = (max(thick) > max(thin_target) > 2.0 < 10.)
     plot, X, Thick, yrange=[0,ymax], $
        xtitle='Distance X from center line (m)', $
        ytitle='Dose (uSv/h)', title="Back Wall Shielding", /nodata
     oplot, X, Thick, color=red
     oplot, X, Thin_target, line=1, color=green
     oplot, X, Thin_target*0.+1.25, line=2, color=yellow
     x0=.4
     y0=.9
     xyouts, x0-.05, y0, /norm, 'Pb shielding thickness = '+$
                      string(Pb, format='(i2)') + 'mm'
     y0 = y0-.05
     xyouts, x0-.05, y0, /norm, $
                        'Distance Z from scatterer to back wall = '+$
                         string(zdist, format='(f5.1)') + 'm'
     y0 = y0-.05
     plots, [x0-.05,x0-.01], [y0,y0], /norm, color=red
     xyouts, x0, y0, /norm, 'Thick scatterer'
     y0 = y0-.05
     plots, [x0-.05,x0-.01], [y0,y0], linestyle=1, /norm, color=green
     xyouts, x0, y0, /norm, 'Thin scatterer'
     y0 = y0-.05
     plots, [x0-.05,x0-.01], [y0,y0], linestyle=2, /norm, color=yellow
     xyouts, x0, y0, /norm, 'Maximum allowed dose'
     openw, 1, file+".txt"
     printf, 1, 'Back wall shielding calculations:'
     printf, 1, 'Pb shielding thickness = '+ string(Pb, format='(i2)') + 'mm'
     printf, 1, 'Distance Z from scatterer to back wall = '+$
                      string(1000*zdist, format='(i4)') + 'mm'
     printf, 1, 'Dose units are uSv/hr'
     printf, 1, '    X Distance  Thick target Thin target'
     printf, 1, '       (m)         dose         dose'
     for j=0, n_elements(X)-1 do begin
        printf, 1, X(j), thick(j), thin_target(j), format='(3f12.3)'
     endfor
     close, 1
  endif else begin
     Z = zend*findgen(100)/99.  ; Distance downstream of scatterer
     brem_side, Z, xdist, Pb, Thick, Thin_target
     ymax = (max(thick) > max(thin_target) > 2.0 < 10.)
     plot, Z, Thick, yrange=[0,ymax], $
        xtitle='Distance Z downstream from scatterer (m)', $
        ytitle='Dose (uSv/h)', title="Side Wall Shielding", /nodata
     oplot, Z, Thick, color=red
     oplot, Z, Thin_target, line=1, color=green
     oplot, Z, Thin_target*0.+1.25, line=2, color=yellow
     x0=.4
     y0=.9
     xyouts, x0-.05, y0, /norm, 'Pb shielding thickness = '+$
                      string(Pb, format='(i2)') + 'mm'
     y0 = y0-.05
     xyouts, x0-.05, y0, /norm, 'Distance X from scatterer to side wall= '+$
                      string(1000*xdist, format='(i4)') + 'mm'
     y0 = y0-.05
     plots, [x0-.05,x0-.01], [y0,y0], /norm, color=red
     xyouts, x0, y0, /norm, 'Thick scatterer'
     y0 = y0-.05
     plots, [x0-.05,x0-.01], [y0,y0], linestyle=1, /norm, color=green
     xyouts, x0, y0, /norm, 'Thin scatterer'
     y0 = y0-.05
     plots, [x0-.05,x0-.01], [y0,y0], linestyle=2, /norm, color=yellow
     xyouts, x0, y0, /norm, 'Maximum allowed dose'
     openw, 1, file+".txt"
     printf, 1, 'Side wall shielding calculations:'
     printf, 1, 'Pb shielding thickness = '+ string(Pb, format='(i2)') + 'mm'
     printf, 1, 'Distance X from beam to inside of shielding = '+$
                      string(1000*xdist, format='(i4)') + 'mm'
     printf, 1, 'Dose units are uSv/hr'
     printf, 1, '    Z Distance  Thick target Thin target'
     printf, 1, '        (m)        dose         dose'
     for j=0, n_elements(Z)-1 do begin
        printf, 1, Z(j), thick(j), thin_target(j), format='(3f12.3)'
     endfor
     close, 1
  endelse
  if (!d.name eq 'Z') then begin
     image = tvrd()
     tvlct, /get, red_lut, green_lut, blue_lut
     write_gif, file+".gif", image, red_lut, green_lut, blue_lut
  endif
endfor
end
