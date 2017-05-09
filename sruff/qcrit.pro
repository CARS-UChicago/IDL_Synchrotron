Function Qcrit, en, elements = elar, dens = ro, dfac = dfac, weights = warr, $
	enuns = dun, formula = form, full = ful, complex = cmp, _extra = _e

;+
; NAME:
;		QCRIT
; VERSION:
;		7.09
; PURPOSE:
;		Calculates Qc for a given material.
; CATEGORY:
;		X-ray calculations.
; CALLING SEQUENCE:
;		Result = QCRIT ([EN,] ELEMENT = ELAR [, optional keywords])
; INPUTS:
;	EN
;		Energy, assumed in the units specified by the variable ENUN in the
;		common block SXR_STUFF, unless specified otherwise by the keyword ENUNS.
;		Needed only if the keywords /FULL and/or /COMPLEX are set.
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
;	/FORMULA
;		Switch.  If set and WEIGHTS (see above) are provided, they're taken as
;		"formula weights".
;	/FULL
;		Switch.  If set, the dispersion correction to Qc is calculated.  In
;		this case, the energy needs to be provided.
;	/COMPLEX
;		Switch.  If set, the complex form of Qc, including the absorption
;		correction, is calculated.  In this case, also, the energy needs to be
;		provided.
;	_EXTRA
;		A formal keyword used to pass keywords to the function DISPER_FUN.  Not
;		to be used directly.  See comment in RESTRICTIONS.
; OUTPUTS:
;		Returns the value(s) of Qc.  If EN is provided the format of the output
;		is same as the format of EN, else it is an energy independent scalar.
;
;		Note:	If /COMPLEX is set, the result is complex, else it is float.
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
;		SXR_STUFF).  Uses calls to ABS_COEFFS, ECONV, ELECOMP, DISPER_FUN and
;		LOAD_ABS_COEFFS from SRUFF_LIB.  Also calls DEFAULT, FPU_FIX, ISNUM and
;		TOLER from MIDL.
; MODIFICATION HISTORY:
;		Created 25-JUL-2002 by Mati Meron as a straightforward modification of
;		the routine DIELECT.
;		Modified 10-JUN-2003 by Mati Meron.  Added keyword DFAC.
;		Modified 15-JUN-2009 by Mati Meron.  Added keyword COMPLEX, for full
;		(complex) calculation.
;-

	common sxr_stuff, sorlist, abctab, enun, denun, curbeam

	on_error, 1
	Load_abs_coeffs
	coc = 8.5300705e-4
	doc = 1.5128681e-5
	ioc = 2.3764076e-005

	enun = Default(dun,enun,/strict)
	dfac = Default(dfac,1.,/dtyp)
	efl = Isnum(en)
	if efl then wen = Econv(en, from = enun, to = denun)

	eli = Elecomp(elar, number = nel)
	if nel gt 1 then begin
		if n_elements(ro) eq 1 then begin
			if n_elements(warr) eq nel then begin
				if keyword_set(form) then twarr = warr*abctab[eli].a $
				else twarr = warr
				tro = ro/total(twarr)*twarr
			endif else message, 'Bad or missing composition data!
		endif else message, 'Bad or missing density!'
	endif else tro = Default(ro,dfac*abctab[eli].ro)
	res = replicate(total(abctab[eli].z*tro/abctab[eli].a),n_elements(en)>1)

	if keyword_set(ful) then begin
		if efl then begin
			for j = 0l, nel - 1 do begin
				k = eli[j]
				ned = abctab[k].edlen
				jmp = abctab[k].edtab[0:ned,0]*abctab[k].edtab[0:ned,1]
				for l = 0l, ned do begin
					q = abctab[k].edtab[l,2]
					t = abctab[k].edtab[l,0]/wen
					if l eq 0 then $
					pres = jmp[0]*t^(1+q)*Disper_fun(0,q,_extra=_e) $
					else pres = pres + jmp[l]*t^(1+q)*Disper_fun(t,q,_extra=_e)
				endfor
				res = (res + doc*tro[j]*pres) > res*Toler(res)
			endfor
		endif else message, "Energy needed for full evaluation!'
	endif

	if keyword_set(cmp) then begin
		if efl then begin
			abc = Abs_coeff(en,elem=elar,dens=ro,dfac=dfac,wei=warr,form=form)
			res = complex(res,-ioc*wen*abc)
		endif else message, "Energy needed for complex part!'
	endif

	if n_elements(res) eq 1 then res = res[0]
	return, FPU_fix(sqrt(coc*res))
end