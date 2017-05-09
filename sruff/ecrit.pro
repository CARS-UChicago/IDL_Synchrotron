Function Ecrit, tet, elements = elar, dens = ro, dfac = dfac, weights = warr, $
	enuns = dun, degrees = deg, formula = form, _extra = _e

;+
; NAME:
;		ECRIT
; VERSION:
;		4.3
; PURPOSE:
;		Calculates Ec for a given material.
; CATEGORY:
;		X-ray calculations.
; CALLING SEQUENCE:
;		Result = ECRIT (TET, [EN,] ELEMENT = ELAR [, optional keywords])
; INPUTS:
;	TET
;		Angle of incidence, in radians, unless /DEGREES is set.
; OPTIONAL INPUT PARAMETERS:
;		None.
; KEYWORD PARAMETERS:
;	ELEMENTS
;		Accepts a list of one or more elements which comprise the target.  Can
;		be provided as character array (each element represented by its
;		chemical symbol) or numeric array (each element represented by its Z).
;	DENS
;		Target density, scalar value.  If provided, overrides the value from
;		the ABCTAB table.  Must be provided if the target is not a pure element
;		(i.e. if ELAR is a vector).
;	DFAC
;		Multiplier for elemental densities, used to account for the possibility
;		that layer densities differ from bulk densities.  Applies only to
;		single element targets.  Ignored if DENS is provided.
;	WEIGHTS
;		In case the element list contains more than one element, the relative
;		partial weights of the elements in the target must be provided as a
;		numeric array, using this keyword.
;	ENUNS
;		Character value, specifies energy units.  Acceptable units are:
;		KEV, EV, NANOMETERS, ANGSTREM.  Only first 2 letters matter.  Default
;		units are specified by the variable ENUN in the common block SXR_STUFF,
;		initially set to the value of DENUN ('keV').  If provided, the value in
;		ENUNS replaces the previous ENUN.
;	/DEGREES
;		Switch.  Specifies that the TET is given in degrees. Default is radians.
;	/FORMULA
;		Switch.  If set and WEIGHTS (see above) are provided, they're taken as
;		"formula weights".
;	_EXTRA
;		A formal keyword used to pass keywords to the function DISPER_FUN.  Not
;		to be used directly.  See comment in RESTRICTIONS.
; OUTPUTS:
;		Returns the value(s) of Ec.  Output format is same as TET.
; OPTIONAL OUTPUT PARAMETERS:
;		None.
; COMMON BLOCKS:
;		SXR_STUFF.  See LOAD_ABS_COEFFS for more information.
; SIDE EFFECTS:
;		If the keyword ENUNS is used, the value of ENUN in the common block
;		SXR_STUFF will be replaced by the new value.
; RESTRICTIONS:
;		None other that the elements specified must be among those in ABCTAB.
; PROCEDURE:
;		Standard evaluation using QCRIT.  See details there.  Calls QCRIT from
;		SRUFF_LIB.  Also calls TOLER and FPU_FIX from MIDL
; MODIFICATION HISTORY:
;		Created 5-AUG-2003 by Mati Meron as a front end to QCRIT.
;-

	common sxr_stuff, sorlist, abctab, enun, denun, curbeam

	on_error, 1
	conv = 12.398424

	if keyword_set(deg) then amult = !dtor else amult = 1.
	qc = (Qcrit(elements = elar, dens = ro, dfac = dfac, weights = warr, $
		enuns = dun, formula = form, _extra = _e))[0]
	res = conv/(4*!pi*sin(amult*tet))*qc

	dif = max(res)
	repeat begin
		pres = res
		pdif = dif
		qc = Qcrit(res,elements = elar, dens = ro, dfac = dfac, $
		weights = warr, enuns = dun, formula = form, /full, _extra = _e)
		res = conv/(4*!pi*sin(amult*tet))*qc
		dif = max(abs(pres - res))
	endrep until dif le 2*Toler() or dif ge pdif

	return, FPU_fix(res)
end