Pro Intprop, r, t, zi = zi, chi = chi, psi = psi, reverse= rev, normalize= norm

;+
; NAME:
;		INTPROP
; VERSION:
;		5.5
; PURPOSE:
;		Calculates refraction in an interface.
; CATEGORY:
;		X-ray specific.
; CALLING SEQUENCE:
;		INTPROP, R [,T], ZI= ZI [CHI=CHI, PSI=PSI, REVERSE=REV, NORMALIZE=NORM]
; INPUTS:
;	R
;		A vector (scalar allowed) of back-propagating (reflected) beam's
;		amplitude.  Complex, in general.
; OPTIONAL INPUT PARAMETERS:
;	T
;		Same as R, for forward propagating amplitude.  If given, must be same
;		length as R.  Defaults to vector of 1s.
; KEYWORD PARAMETERS:
;	ZI
;		Vector of Zi values, see Refraction in Layered Matter (RLM).  Same
;		length as R.
;	CHI
;		A vector of Chi values (see RLM).  Optional, if given same length as R.
;	PSI
;		A vector of Psi values (see RLM).  Optional, if given same length as R.
;	/REVERSE
;		Switch.  If set, backwards propagation is performed.
;	/NORMALIZE
;		Switch.  If set, the result is normalized to value(s) of 1 for T.
; OUTPUTS:
;		On return, R and T (if given) contain the propagated values.
; OPTIONAL OUTPUT PARAMETERS:
;		None.
; COMMON BLOCKS:
;		None.
; SIDE EFFECTS:
;		R and T change during the calculation.
; RESTRICTIONS:
;		All vectors must be of same length.
; PROCEDURE:
;		Based on RLM.  Calls CALCTYPE and DEFAULT from MIDL.
; MODIFICATION HISTORY:
;		Created 5-JUN-2006 by Mati Meron.
;-

	on_error, 1

	len = n_elements(r)
	if len eq 0 or len ne n_elements(zi) then message, 'Bad data!'
	typ = Calctype(r,zi)
	def = make_array(size=size(r),type=typ,val=1)
	t = Default(t,def)
	chi = Default(chi,def)
	psi = Default(psi,def)
	if n_elements(chi) ne len or n_elements(psi) ne len $
	or n_elements(t) ne len then message,'Bad extra data!'

	if keyword_set(rev) then begin
		denom = psi*(1 + zi)
		rtem = (chi*psi*r + chi*zi*t)/denom
		t = (psi*zi*r + t)/denom
		r = rtem
	endif else begin
		denom = chi*(1 - zi)
		rtem = (r - chi*zi*t)/denom
		t = (-psi*zi*r + chi*psi*t)/denom
		r = rtem
	endelse

	if keyword_set(norm) then begin
		r = r/t
		t = def
	endif

	return
end