; $Id$

pro brem_side, z, x, pb, thick, thin

; NAME:
;	BREM_SIDE
;
; PURPOSE:
;       This procedure computes the dose from secondary bremsstrahlung thru 
;       the side wall of an enclosure or beam transport. 
;
; CATEGORY:
;	Synchrotron beamline shielding
;
; CALLING SEQUENCE:
;	ROUTINE_NAME, Z, X, Pb, Thick, Thin
;
;
; INPUTS:
;	Z:      The distance downstream from the scatterer at which to
;               calculate the dose. This is the distance measured along the
;               direction of beam propogation, and is not the total distance
;               to the observation point.  Z can be an array, in which case the
;               output values (Thick and Thin) will also be arrays. 
;               Units=meters
;
;	X:      The lateral distance from the beam to the inside of the
;               shielding surface.  The observation point is an additional
;               0.3m outside the shielding surface.
;               Units = meters.  X can be an array, in which case it must have
;               the same dimensions as Z. If X is a scalar then X is
;               the same for all values of Z.
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
;       Calculate the dose from 0 to 5 meters downstream of a scatterer
;       through a 11mm Pb beampipe wall located 5cm from the direct beam.
;       Z = findgen(50)/10.
;	BREM_SIDE, Z, .05, 11., Thick, Thin
;       PLOT, Z, Thick
;       OPLOT, Z, Thin
;
; MODIFICATION HISTORY:
; 	Written by:	Mark Rivers, based heavily upon Gerd Rosenbaum's
;	June, 1995      program.
;-

; Check dimensions of X and Z are consistent
if ((n_elements(X) gt 1) and (n_elements(X) ne n_elements(Z))) then $
    message, 'Dimensions of X and Z must agree if X is an array'

m1 = Pb     ; lead thickness (mm)
m2 = X      ; lateral wall distance (m)
m3 = 0.0484*m1*(1. - 0.254*exp(-m1/20.4))
m4 = 12.7 * (1.-exp(-m1/21.5))

c2 = Z
c0 = atan(m2+0.3,c2)/!dtor
c1 = (m2+0.3)/sin(c0*!dtor)
c4 = c2
; Compute thick target value
for i=0, n_elements(c2)-1 do begin
   if (c1(i) lt 0.1) then c4(i)=1. $
   else if (c1(i) lt 1.) then c4(i)=0.1/c1(i) $
   else c4(i)=0.1/c1(i)^2
endfor
c5 = 247.*exp(-c0/38.8)*c4*exp(-m3*(1.-exp(-c0/m4))/sin(c0*!dtor))
thick = c5
; Compute thin target value
for i=0, n_elements(c2)-1 do begin
   if (c1(i) lt 0.1) then c4(i)=1. $
   else if (c1(i) lt .3) then c4(i)=0.1/c1(i) $
   else c4(i)=0.027/c1(i)^2
endfor
c5 = 10^(4.42 - 1.20*alog10(c0)-.221*(alog10(c0))^2-0.124*(alog10(c0))^3)*c4 $
     *exp(-m3*(1-exp(-c0/m4))/sin(c0*!dtor))
thin = c5
end
