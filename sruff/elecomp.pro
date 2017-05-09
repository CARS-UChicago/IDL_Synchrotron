Function Elecomp, elar, number = nel

;+
; NAME:
;		ELECOMP
; VERSION:
;		4.2
; PURPOSE:
;		Translates a list of elements into a list of their respective indexes
;		in the table ABCTAB (common block SXR_STUFF).
; CATEGORY:
;		X-ray calculations.
; CALLING SEQUENCE:
;		Result = ELECOMP (ELAR [, NUMBER = NEL])
; INPUTS:
;	ELAR
;		List of elements, provided either as character array (with each element
;		represented by its chemical symbol) or as a numeric array (with each
;		element represented by its Z value).
; OPTIONAL INPUT PARAMETERS:
;		None.
; KEYWORD PARAMETERS:
;	NUMBER
;		Optional output parameters, see below.
; OUTPUTS:
;		Returns an integer array, where each entry is the index of the
;		corresponding element from ELAR in the table ABCTAB.
; OPTIONAL OUTPUT PARAMETERS:
;	NUMBER
;		The name of a variable to receive the number of elements in ELAR.
; COMMON BLOCKS:
;		SXR_STUFF.  Contains (second entry) the table ABCTAB.
; SIDE EFFECTS:
;		None.
; RESTRICTIONS:
;		None other that the elements specified must be among those in ABCTAB.
; PROCEDURE:
;		Straightforward.  Using TYPE and STRMATCH_MM from MIDL.
; MODIFICATION HISTORY:
;		Created 1-MARCH-1993 by Mati Meron.
;		Modified 15-SEP-2001 by Mati Meron.  Verified WINDOWS compatibility.
;-

	common sxr_stuff, sorlist, abctab, enun, denun, curbeam

	on_error, 1
	Load_abs_coeffs

	nel = n_elements(elar)
	elin = intarr(nel)
	if Type(elar) eq 7 then begin
		for j = 0, nel - 1 do elin[j] = StrMatch_mm(elar[j],abctab.csym,2)
	endif else for j = 0, nel - 1 do elin[j] = where(abctab.z eq elar[j])
	ver = where(elin eq -1, vnum)
	if vnum gt 0 then message, 'Some data is missing!'

	return, elin
end