function mono_crystal, $
         angle=angle, $
         energy=energy, $
         silicon=silicon, $
         germanium=germanium, $
         hkl=hkl

; This function computes the angle or energy of a monochromator crystal
; One of the keyword parameters "energy" or "angle" must be specified.  
; If energy is specified then the function returns the angle, and if angle is
; specified then the function returns the energy.
; Units:
;   Angle: degreess
;   Energy: keV
; The silicon and germanium keywords are mutually exclusive.  
;   The default is silicon.
; The default value of hkl is [1,1,1]
; Examples:
; To compute the angle of Si(111) at 15 keV
;   print, mono_crystal(energy=15)
; To compute the energy of Ge(220) at 23 degrees
;   print, mono_crystal(angle=23, /germanium, hkl=[2,2,0])
asi = 5.43094
age = 5.32

if n_elements(germanium) ne 0 then a0 = age else a0 = asi
if n_elements(hkl) eq 0 then hkl=[1,1,1]
shkl = sqrt(hkl(0)^2 + hkl(1)^2 + hkl(2)^2)
if n_elements(energy) ne 0 then begin
   lambda = 12.398/energy
   theta = asin(lambda * shkl/(2. * a0))
   return, theta/!dtor
endif else begin
   lambda = sin(angle*!dtor) * 2. * a0 / shkl
   return, 12.398/lambda
endelse
end
