Function Dielect, en, elements= elar, dens= ro, dfac= dfac, weights= warr, $
	enuns = dun, formula = form, _extra = _e

;+
; NAME:
;		DIELECT
; VERSION:
;		4.2
; PURPOSE:
;		Calculates the energy dependent part of the dielectric constant at
;		X-ray energies, i.e. 1 - eps(e).
; CATEGORY:
;		X-ray calculations.
; CALLING SEQUENCE:
;		Result = DIELECT (EN, ELEMENT = ELAR [, optional keywords])
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
;	/FORMULA
;		Switch.  If set and WEIGHTS (see above) are provided, they're taken as
;		"formula weights".
;	_EXTRA
;		A formal keyword used to pass keywords to the function DISPER_FUN.  Not
;		to be used directly.  See comment in RESTRICTIONS.
; OUTPUTS:
;		Returns 1 - (dielectric coefficient) for all the energies in EN.
;		Output form is same as EN.
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
;		Standard evaluation using densities from ABCTAB (in common block
;		SXR_STUFF) and absorption coefficients provided by ABS_COEFF.  Uses
;		calls to ECONV, ELECOMP, DISPER_FUN and LOAD_ABS_COEFFS.  Also calls
;		DEFAULT, FPU_FIX and TOLER from MIDL.
; MODIFICATION HISTORY:
;		Created 1-MARCH-1993 by Mati Meron.
;		Modified 1-JUNE-1994 by Mati Meron.  Added the option of user specified
;		density.
;		Modified 15-SEP-2001 by Mati Meron.  Automated the call to
;		LOAD_ABS_COEFFS and verified WINDOWS compatibility.
;		Modified 25-JUL-2002 by Mati Meron.  Added keyword FORMULA.
;		Modified 10-JUN-2003 by Mati Meron.  Added keyword DFAC.
;-

	common sxr_stuff, sorlist, abctab, enun, denun, curbeam

	on_error, 1
	Load_abs_coeffs
	rpc = 8.3035847e-4
	dpc = 1.2562227e-8
	ipc = 1.9732700e-8

	dfac = Default(dfac,1.,/dtyp)
	enun = Default(dun,enun,/strict)
	in_units = enun
	out_units = denun
	ien = 1./Econv(en, from = in_units, to = out_units)
	ipart = ipc*Abs_coeff(en,elem=elar,dens=ro,dfac=dfac,wei=warr,form=form)

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

	rpart = rpc*total(abctab[eli].z*tro/abctab[eli].a)
	for j = 0l, nel - 1 do begin
		k = eli[j]
		ned = abctab[k].edlen
		jmp = abctab[k].edtab[0:ned,0]*abctab[k].edtab[0:ned,1]
		for l = 0l, ned do begin
			q = abctab[k].edtab[l,2]
			t = abctab[k].edtab[l,0]*ien
			if l eq 0 then pres = jmp[0]*t^(1+q)*Disper_fun(0,q,_extra=_e) $
			else pres = pres + jmp[l]*t^(1+q)*Disper_fun(t,q,_extra=_e)
		endfor
		rpart = (rpart + dpc*tro[j]*pres) > rpart*Toler(rpart)
	endfor

	return, FPU_fix(ien*complex(rpart*ien, -ipart))
end