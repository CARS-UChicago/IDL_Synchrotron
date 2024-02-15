Function Abs_coeff, en, elements= elar, dens= ro, dfac= dfac, weights= warr, $
     enuns = dun, mass = mas, formula = form

;+
; NAME:
;		ABS_COEFF
; VERSION:
;		4.3
; PURPOSE:
;		Calculates absorption coefficients (NOT mass absorption) for a given
;		target and X-ray energy.
; CATEGORY:
;		X-ray calculations.
; CALLING SEQUENCE:
;		Result = ABS_COEFF (EN,  ELEMENT = ELAR [, optional keywords])
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
;		Density multiplier, used to account for the possibility that layer
;		densities differ from bulk densities.  Default value is 1.
;	WEIGHTS
;		In case the element list contains more than one element, the relative
;		partial weights of the elements in the target must be provided as a
;		numeric array, using this keyword.
;	ENUNS
;		Character value, specifies energy units.  Acceptable units are:
;		KEV, EV, NANOMETERS, ANGSTREM.  Only first 2 letters matter.  Default
;		units are specified by the variable ENUN in the common block SXR_STUFF,
;		initially set to the value of DEFUN ('keV').  If provided, the value in
;		ENUNS replaces the previous ENUN.
;	/MASS
;		Switch.  If set, mass absorption coefficient is returned.
;	/FORMULA
;		Switch.  If set and WEIGHTS (see above) are provided, they're taken as
;		"formula weights".
; OUTPUTS:
;		Returns the target's absorption coefficient, for each energy in EN.
;		Format of result is same as the format of EV.  The units are 1/cm,
;		unless MASS is set, in which case the the units are cm^2/gr.
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
;		Spline evaluation using the spline coefficients in the table ABCTAB
;		(common block SXR_STUFF).  Uses calls to ECONV, ELECOMP and, if needed,
;		LOAD_ABS_COEFFS.  Also calls DEFAULT, FPU_FIX and SPLIN_EVAL from MIDL.
; MODIFICATION HISTORY:
;		Created 1-MAR-1993 by Mati Meron.
;		Modified 1-JUN-1994 by Mati Meron.  Added the option of user specified
;		density and the keyword MASS.
;		Modified 15-SEP-2001 by Mati Meron.  Automated the call to
;		LOAD_ABS_COEFFS and verified WINDOWS compatibility.
;		Modified 25-JUL-2002 by Mati Meron.  Added keyword FORMULA.
;		Modified 10-JUN-2003 by Mati Meron.  Added keyword DFAC.
;-

	common sxr_stuff, sorlist, abctab, enun, denun, curbeam

	on_error, 1
	Load_abs_coeffs
	dfac = Default(dfac,1.,/dtyp)
	enun = Default(dun,enun,/strict)
	rmax = alog(1e6)
	in_units = enun
	out_units = denun

	logen = alog(Econv(en, from = in_units, to = out_units))
	eli = Elecomp(elar, number = nel)
	if nel gt 1 then begin
		if n_elements(ro) eq 1 then begin
			if n_elements(warr) eq nel then begin
				if keyword_set(form) then twarr = warr*abctab[eli].a $
				else twarr = warr
				tro = dfac*ro/total(twarr)*twarr
			endif else message, 'Bad or missing composition data!
		endif else message, 'Bad or missing density!'
	endif else tro = dfac*Default(ro,abctab[eli].ro)
	if keyword_set(mas) then tro = tro/total(tro)

	res = 0.*logen
	for j = 0, nel - 1 do begin
		k = eli[j]
		res = res + $
		tro[j]*exp(Splin_eval(logen,abctab[k].cotab[0:abctab[k].colen,*]) <rmax)
	endfor

	return, FPU_fix(res)
end