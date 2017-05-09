Pro Homprop, r, t, kni = kni, dist = dist, reverse= rev, normalize= norm

;+
; NAME:
;		HOMPROP
; VERSION:
;		5.5
; PURPOSE:
;		Calculates propagation in homogenous matter..
; CATEGORY:
;		X-ray specific.
; CALLING SEQUENCE:
;		HOMPROP, R [,T], KNI= KNI [REVERSE=REV, NORMALIZE=NORM]
; INPUTS:
;	R
;		A vector (scalar allowed) of back-propagating (reflected) beam's
;		amplitude.  Complex, in general.
; OPTIONAL INPUT PARAMETERS:
;	T
;		Same as R, for forward propagating amplitude.  If given, must be same
;		length as R.  Defaults to vector of 1s.
; KEYWORD PARAMETERS:
;	KNI
;		Vector of k*ni values, see Refraction in Layered Matter (RLM).  Same
;		length as R.
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
	if len eq 0 or len ne n_elements(kni) then message, 'Bad data!'
	typ = Calctype(0.,r,kni)
	def = make_array(size=size(r),type=typ,val=1)
	t = Default(t,def)
	if n_elements(t) ne len then message,'Bad extra data!'
	if keyword_set(rev) then z = -dist[0] else z = dist[0]
	eiphi = exp(complex(0,z)*kni)

	r = r/eiphi
	t = t*eiphi

	if keyword_set(norm) then begin
		r = r/t
		t = def
	endif

	return
end