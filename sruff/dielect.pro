Function Dielect, en, elements = elarr, dens = ro, weights = warr, enuns = dun

;+
; NAME:
;	DIELECT
; PURPOSE:
;	Calculates the energy dependent part of the dielectric constant at 
;	X-ray energies, i.e. 1 - eps(e).
; CATEGORY:
;	X-ray calculations.
; CALLING SEQUENCE:
;	Result = DIELECT (EN, ELEMENT = ELARR [, optional keywords])
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
;	Character value, specifies energy units.  Acceptable units are:
;	KEV, EV, NANOMETERS, ANGSTREM.  Only first 2 letters matter.  Default
;	units are specified by the variable ENUN in the common block SXR_STUFF,
;	initially set to the value of DEFUN ('keV').  If provided, the value in
;	ENUNS replaces the previous ENUN.
; OUTPUTS:
;	Returns 1 - (dielectric coefficient) for all the energies in EN.  
;	Output form is same as EN.
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
;	Standard evaluation using densities from ABCTAB (in common block 
;	SXR_STUFF) and absorption coefficients provided by ABS_COEFF.  Uses 
;	calls to ECONV, ELECOMP and DISPER_FUN.
; MODIFICATION HISTORY:
;	Created 1-MARCH-1993 by Mati Meron.
;	Modified 1-JUNE-1994 by Mati Meron.  Added the option of user specified
;	density.
;-

    common sxr_stuff, sorlist, abctab, enun, denun, curbeam

    load_abs_coeffs
    on_error, 1
    if n_elements(abctab) eq 0 then message, 'Tables not loaded!'
    enun = Default(dun,enun,/strict)
    rpc = 8.3035847e-4
    dpc = 1.2562227e-8
    ipc = 1.9732700e-8

    ien = 1./Econv(en, from = enun, to = denun)
    ipart = ipc*Abs_coeff(en, elements = elarr, dens = ro, weights = warr)

    eli = Elecomp(elarr, number = nel)
    if nel eq 1 then tro = Default(ro,abctab(eli).ro) else $
    if n_elements(ro) ne 1 then message, 'Bad or missing density!' else $
    if n_elements(warr) ne nel then message, 'Missing composition data!' else $
    tro = ro/total(warr)*warr

    rpart = rpc*total(abctab(eli).z*tro/abctab(eli).a)
    for j = 0l, nel - 1 do begin
	k = eli(j)
	ned = abctab(k).edlen
	jmp = abctab(k).edtab(0:ned,0)*abctab(k).edtab(0:ned,1)
	for l = 0l, ned do begin
	    q = abctab(k).edtab(l,2)
	    t = abctab(k).edtab(l,0)*ien
	    if l eq 0 then pres = jmp(0)*t^(1+q)*Disper_fun(0,q) else $
	    pres = pres + jmp(l)*t^(1+q)*Disper_fun(t,q)
	endfor
	rpart = rpart + dpc*tro(j)*pres
    endfor

    return, ien*complex(rpart*ien, -ipart)
end
