Function Abs_length, en, elements= elar, dens= ro, dfac= dfac, weights= warr, $
	formula= form

;+
; NAME:
;		ABS_LENGTH
; VERSION:
;		4.3
; PURPOSE:
;		Calculates absorption lengths for a given target and X-ray energy.
; CATEGORY:
;		X-ray calculations.
; CALLING SEQUENCE:
;		Result = ABS_LENGTH (EN,  ELEMENT = ELAR [, optional keywords])
; INPUTS:
;	EN
;		Energy, assumed in the units specified by the variable ENUN in the
;		common block SXR_STUFF, unless specified otherwise by the keyword ENUNS.
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
;	/FORMULA
;		Switch.  If set and WEIGHTS (see above) are provided, they're taken as
;		"formula weights".
; OUTPUTS:
;		Returns the target's absorption length, for each energy in EN.
;		Format of result is same as the format of EN.  The units are cm,
; OPTIONAL OUTPUT PARAMETERS:
;		None.
; COMMON BLOCKS:
;		None.
; SIDE EFFECTS:
;		None.
; RESTRICTIONS:
;		None other then those specified for ABS_COEFFS.
; PROCEDURE:
;		Front end for ABS_COEFFS, in SRUFF.  Also calls FPU_FIX from MIDL.
; MODIFICATION HISTORY:
;		Created 25-JUL-2002 by Mati Meron as a straightforward application of
;		ABS_COEFFS.
;		Modified 10-JUN-2003 by Mati Meron.  Added keyword DFAC.
;-

	return, $
	FPU_fix(1./Abs_coeff(en,elem=elar,dens=ro,dfac=dfac,wei=warr,form=form))
end