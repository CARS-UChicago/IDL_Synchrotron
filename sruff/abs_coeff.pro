Function Abs_coeff, en, elements = elarr, dens = ro, weights = warr, $
    enuns = dun, mass = mas

;+
; NAME:
;	ABS_COEFF
; PURPOSE:
;	Calculates absorption coefficients (NOT mass absorption) for a given
;	target and X-ray energy.
; CATEGORY:
;	X-ray calculations.
; CALLING SEQUENCE:
;	Result = ABS_COEFF (EN,  ELEMENT = ELARR [, optional keywords])
; INPUTS:
;    EN
;	Energy, assumed in the units specified by the variable ENUN in the 
;	common block SXR_STUFF, unless specified otherwise by the keyword ENUNS.
; OPTIONAL INPUT PARAMETERS:
;	None.
; KEYWORD PARAMETERS:
;    ELEMENTS
;	Accepts a list of one or more elements which comprise the target.  Can
;	be provided as character array (each element represented by its 
;	chemical symbol) or numeric array (each element represented by its Z).
;    DENS
;	Target density, scalar value.  If provided, overrides the value from 
;	the ABCTAB table.  Must be provided if the target is not a pure element
;	(i.e. if ELARR is a vector).
;    WEIGHTS
;	In case the element list contains more than one element, the relative 
;	partial weights of the elements in the target must be provided as a
;	numeric array, using this keyword.
;    ENUNS
;       Character value, specifies energy units.  Acceptable units are:
;       KEV, EV, NANOMETERS, ANGSTREM.  Only first 2 letters matter.  Default
;	units are specified by the variable ENUN in the common block SXR_STUFF,
;	initially set to the value of DEFUN ('keV').  If provided, the value in
;	ENUNS replaces the previous ENUN.
;    /MASS
;	Switch.  If set, mass absorption coefficient is returned.
; OUTPUTS:
;	Returns the target's absorption coefficient, for each energy in EN.
;	Format of result is same as the format of EV.  The units are 1/cm, 
;	unless MASS is set, in which case the the units are cm^2/gr.
; OPTIONAL OUTPUT PARAMETERS:
;	None.
; COMMON BLOCKS:
;	SXR_STUFF.  See LOAD_ABS_COEFFS for more information.
; SIDE EFFECTS:
;	If necessary the common block SXR_STUFF will be initialized via a call 
;       to LOAD_ABS_COEFFS.
;	If the keyword ENUNS is used, the value of ENUN in the common block 
;	SXR_STUFF will be replaced by the new value.
; RESTRICTIONS:
;	The environment variable IDL_ABS_COEFFS should point to the XDR binary
;       file of absorption coefficients.  This file is read by LOAD_ABS_COEFFS.
; PROCEDURE:
;	Spline evaluation using the spline coefficients in the table ABCTAB 
;	(common block SXR_STUFF).  Uses calls to ECONV and ELECOMP.  Also calls 
;	DEFAULT and SPLIN_EVAL from MIDL.
; MODIFICATION HISTORY:
;	Created 1-MARCH-1993 by Mati Meron.
;	Modified 1-JUNE-1994 by Mati Meron.  Added the option of user specified
;	density and the keyword MASS.
;-

    common sxr_stuff, sorlist, abctab, enun, denun, curbeam

    load_abs_coeffs
    on_error, 1
    if n_elements(abctab) eq 0 then message, 'Tables not loaded!'
    enun = Default(dun,enun,/strict)
    rmax = alog(1e6)

    logen = alog(Econv(en, from = enun, to = denun))
    eli = Elecomp(elarr, number = nel)
    if nel eq 1 then tro = Default(ro,abctab(eli).ro) else $
    if n_elements(ro) ne 1 then message, 'Bad or missing density!' else $
    if n_elements(warr) ne nel then message, 'Missing composition data!' else $
    tro = ro/total(warr)*warr
    if keyword_set(mas) then tro = tro/total(tro)

    res = 0.*logen
    for j = 0, nel - 1 do begin
	k = eli(j)
	res = res + $
	tro(j)*exp(Splin_eval(logen,abctab(k).cotab(0:abctab(k).colen,*)) <rmax)
    endfor

    return, res
end
