; Program to calculate and plot bremsstrahlung dose thru back walls.

if (!d.name eq 'PS') then begin
  device, xsize=10, ysize=7, /inch, /land, file='brem_back_main.ps'
  !p.font=0
  !p.thick=3
  !p.charsize=1.5
  chsize=1.5
endif else begin
  !p.font=-1
  !p.charsize=1.5
  !p.thick=1
  chsize=1.5
endelse

X = 0.1+findgen(30)/10.  ; Lateral distance to inside of shielding
Z = 1             ; Distance downstream of scatterer
Pb = 50.            ; Pb thickness

brem_back, Z, X, Pb, Thick, Thin

plot, X, Thick, yrange=[0,max(thick) > max(thin)], $
        xtitle='Distance X from center line (m)', $
        ytitle='Dose (uSv/h)'
oplot, X, Thin, line=1
oplot, X, Thin*0.+1.25, line=2
x0=.4
y0=.9
xyouts, x0-.05, y0, /norm, 'Pb shielding thickness = '+$
                      string(Pb, format='(i2)') + 'mm'
y0 = y0-.05
xyouts, x0-.05, y0, /norm, 'Distance Z from scatterer to back wall = '+$
                      string(Z(0), format='(f5.1)') + 'm'
y0 = y0-.05
plots, [x0-.05,x0-.01], [y0,y0], /norm
xyouts, x0, y0, /norm, 'Thick scatterer'
y0 = y0-.05
plots, [x0-.05,x0-.01], [y0,y0], line=1, /norm
xyouts, x0, y0, /norm, 'Thin scatterer'
y0 = y0-.05
plots, [x0-.05,x0-.01], [y0,y0], line=2, /norm
xyouts, x0, y0, /norm, 'Maximum allowed dose'
end
