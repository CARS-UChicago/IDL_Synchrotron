; Program to calculate and plot bremsstrahlung dose thru lateral walls.

if (!d.name eq 'PS') then begin
  device, xsize=10, ysize=7, /inch, /land, file='brem_side.ps'
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

Z = findgen(50)/10.  ; Distance downstream of scatterer
X = .03             ; Lateral distance to inside of shielding
Pb = 19.            ; Pb thickness

brem_side, Z, X, Pb, Thick, Thin

plot, Z, Thick, yrange=[0,max(thick) > max(thin)], $
        xtitle='Distance Z downstream from scatterer (m)', $
        ytitle='Dose (uSv/h)'
oplot, Z, Thin, line=1
oplot, Z, Thin*0.+1.25, line=2
x0=.4
y0=.9
xyouts, x0-.05, y0, /norm, 'Pb shielding thickness = '+$
                      string(Pb, format='(i2)') + 'mm'
y0 = y0-.05
xyouts, x0-.05, y0, /norm, 'Distance X from direct beam to inside of shielding = '+$
                      string(1000*X(0), format='(i4)') + 'mm'
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
