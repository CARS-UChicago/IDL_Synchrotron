Function Thruput, ene, filters = filt, filthicks = flth, $
	mirrors = mirr, mirangles = mran, target = targ, tarthick = tath, _extra= _e

;+
; NAME:
;		THRUPUT
; VERSION:
;		6.0
; PURPOSE:
;		Calculates the throughput of a combination of x-ray mirrors and filters.
; CATEGORY:
;		X-ray calculations.
; CALLING SEQUENCE:
;		Result = THRUPUT (ENE [, FILTERS = FILT, FILTHICKS = FLTH] $
;				[, MIRRORS = MIRR, MIRANGLES = MRAN] $
;				[, TARGET = TARG, TARTHICK = TATH]
; INPUTS:
;	ENE
;		Energy, in keV.  Mandatory.
;	FILTERS														| If given,
;		Character scalar or vector, a list of (single element)	| these two must
;		filters.												| have same
;	FILTHICKS													| number of
;		Numeric scalar or vector, a list of filter thicknesses,	| elements
;		in mm.
;	MIRRORS														| If given,
;		Character scalar or vector, a list of (single element)	| these two must
;		mirror coatings.										| have same
;	MIRANGLES													| number of
;		Numeric scalar or vector, a list mirror angles in mrad.	| elements.
;	TARGET
;		Character scalar the target element.
;	TARTHICK
;		Numeric scalar, the target thickness in mm.  Optional, even if target
;		is given.  Taken as infinite (i.e. full absorption) if not given.
;	_EXTRA
;		A formal keyword used to pass keywords to the function MIRROR.  Not to
;		be used directly.  See comment in RESTRICTIONS.
; OUTPUTS:
;		Returns the fraction of the beam passed through the filters and/or
;		mirrors and (if TARGET and TARTHICK are given) absorbed in the target,
;		for all the energies in ENE.  Output form is same as ENE.
; OPTIONAL OUTPUT PARAMETERS:
;		None.
; COMMON BLOCKS:
;		None.
; SIDE EFFECTS:
;		None.
; RESTRICTIONS:
;		None other than the restrictions of the called functions.
; PROCEDURE:
;		Straightforward.  Calls ABS_COEFF and MIRROR from SRUFF_LIB.  Also
;		calls CALCTYPE, CAST and TOLER from MIDL.
; MODIFICATION HISTORY:
;		Created 5-FEB-2007 by Mati Meron.
;-

	on_error, 1
	amult = 1d-3
	fmult = 1d-1

	eps = Toler()
	typ = Calctype(ene,0.)
	res = 0.*ene + 1.

	nf = n_elements(filt)
	if nf gt 0 then begin
		if n_elements(flth) eq nf then begin
			for i = 0l, nf-1 do $
			res = res*exp(-fmult*flth[i]*Abs_coeff(ene>eps,elem=filt[i]))
		endif else message, 'Inconsistent filter data!'
	endif

	nm = n_elements(mirr)
	if nm gt 0 then begin
		if n_elements(mran) eq nm then begin
			for i = 0l, nm-1 do $
			res = res*Mirror(ene>eps,amult*mran[i],elem=mirr[i],_extra=_e)
		endif else message, 'Inconsistent mirror data!'
	endif

	nt = n_elements(targ)
	if nt eq 1 then begin
		if n_elements(tath) eq 1 then begin
			res = res*(1 - exp(-fmult*tath[0]*Abs_coeff(ene>eps,elem=targ[0])))
		endif else if n_elements(tath) ne 0 then $
		message, 'Only single target thickness allowed!'
	endif else if nt ne 0 then message, 'Only single target allowed'

	return, Cast(res,typ,typ,/fix)
end