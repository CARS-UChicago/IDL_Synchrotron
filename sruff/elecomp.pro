Function Elecomp, elarr, number = nel

;+
; NAME:
;	ELECOMP
; PURPOSE:
;	Translates a list of elements into a list of their respective indexes
;	in the table ABCTAB (common block SXR_STUFF).
; CATEGORY:
;	X-ray calculations.
; CALLING SEQUENCE:
;	Result = ELECOMP (ELARR [, NUMBER = NEL])
; INPUTS:
;    ELARR
;	List of elements, provided either as character array (with each element
;	represented by its chemical symbol) or as a numeric array (with each
;	element represented by its Z value).
; OPTIONAL INPUT PARAMETERS:
;	None.
; KEYWORD PARAMETERS:
;    NUMBER
;	Optional output parameters, see below.
; OUTPUTS:
;	Returns an integer array, where each entry is the index of the 
;	corresponding element from ELARR in the table ABCTAB.
; OPTIONAL OUTPUT PARAMETERS:
;    NUMBER
;	The name of a variable to receive the number of elements in ELARR.
; COMMON BLOCKS:
;	SXR_STUFF.  Contains (second entry) the table ABCTAB.
; SIDE EFFECTS:
;	If necessary the common block SXR_STUFF will be initialized via a call 
;       to LOAD_ABS_COEFFS.
; RESTRICTIONS:
;	The environment variable IDL_ABS_COEFFS should point to the XDR binary
;       file of absorption coefficients.  This file is read by LOAD_ABS_COEFFS.
; PROCEDURE:
;	Straightforward.  Using TYPE and STRMATCH from MIDL.
; MODIFICATION HISTORY:
;	Created 1-MARCH-1993 by Mati Meron.
;-

    common sxr_stuff, sorlist, abctab, enun, denun, curbeam

    load_abs_coeffs
    on_error, 1
    if n_elements(abctab) eq 0 then message, 'Tables not loaded!'

    nel = n_elements(elarr)
    elin = intarr(nel)
    if Type(elarr) eq 7 then begin
	for j = 0, nel - 1 do elin(j) = StrMatch_mm(elarr(j),abctab.csym,2)
    endif else for j = 0, nel - 1 do elin(j) = where(abctab.z eq elarr(j))
    ver = where(elin eq -1, vnum)
    if vnum gt 0 then message, 'Some data is missing!'

    return, elin
end
