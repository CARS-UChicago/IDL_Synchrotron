pro brem_back, z, x, pb, thick, thin

; NAME:
;	BREM_BACK
;
; PURPOSE:
;       This procedure computes the dose from secondary bremsstrahlung thru 
;       the back wall of an enclosure or beam transport. 
;
; CATEGORY:
;	Synchrotron beamline shielding
;
; CALLING SEQUENCE:
;	ROUTINE_NAME, Z, X, Pb, Thick, Thin
;
;
; INPUTS:
;	Z:      The distance from the scatterer to the back wall. The
;               observation point is an additional 0.3 meters behind the 
;               back wall. Z can be an array, in which case the number of
;               elements of Z will be the first dimension of the output arrays.
;               Units=meters
;
;	X:      The lateral distance from the beam to the observation point.
;               Units = meters.  X can be an array, in which case the output
;               values (Thick and Thin) will also be arrays.
;               Units=meters
;
;       Pb:     The thickness of the Pb shielding in mm.
;
; OUTPUTS:
;	Thick:  The dose due to a thick target (5x5x30cm). uSv/hr
;
;	Thin:   The dose due to a thin target (5x5x5cm). uSv/hr
;   
; PROCEDURE:
;	Copy Gerd's documentation here.
;
; EXAMPLE:
;       Calculate the dose from 0 to 3 meters laterally along the back wall
;       from a scatterer located 4 meters upstream of the back wall. The
;       backwall Pb thickness is 19mm.
;       X = FINDGEN(30)/10.
;	BREM_BACK, 5, X, 19., Thick, Thin
;       PLOT, X, Thick
;       OPLOT, X, Thin
;
; MODIFICATION HISTORY:
; 	Written by:	Mark Rivers, based heavily upon Gerd Rosenbaum's
;	June, 1995      program.
;-

; Check dimensions of X and Z are consistent

m1 = Pb     ; lead thickness (mm)
m2 = Z      ; lateral wall distance (m)
m3 = 0.0484*m1*(1. - 0.254*exp(-m1/20.4))
m4 = 12.7 * (1.-exp(-m1/21.5))

c2 = X
c0 = atan(c2,m2+0.3)/!dtor
c1 = (m2+0.3)/cos(c0*!dtor)
c4 = c2
; Compute thick target value
for i=0, n_elements(c2)-1 do begin
   if (c1(i) lt 0.1) then c4(i)=1. $
   else if (c1(i) lt 1.) then c4(i)=0.1/c1(i) $
   else c4(i)=0.1/c1(i)^2
endfor
c5 = 247.*exp(-c0/38.8)*c4*exp(-m3*(1.-exp(-(90.-abs(c0))/m4))/cos(c0*!dtor))
thick = c5
; Compute thin target value
for i=0, n_elements(c2)-1 do begin
   if (c1(i) lt 0.1) then c4(i)=1. $
   else if (c1(i) lt .3) then c4(i)=0.1/c1(i) $
   else c4(i)=0.027/c1(i)^2
endfor
c5 = 10^(4.42 - 1.20*alog10(c0)-.221*(alog10(c0))^2-0.124*(alog10(c0))^3)*c4 $
     *exp(-m3*(1-exp(-(90.-abs(c0))/m4))/cos(c0*!dtor))
thin = c5
end
