Function Mirror, en, theta, elements = elarr, num_elems = numels, $
    dens = rarr, weights = warr, thicks = tarr, enuns = dun, $
    degrees = deg, centimeter = cm, transmition = tran, coherent = coh

;+
; NAME:
;	MIRROR
; PURPOSE:
;	Calculates the reflection (or transmition) coefficient of an X-ray 
;	mirror.
; CATEGORY:
;	X-ray calculations.
; CALLING SEQUENCE:
;	Result = MIRROR (EN, THETA, ELEMENTS = ELARR [, optional keywords])
; INPUTS:
;    EN
;	Energy, assumed in the units specified by the variable ENUN in the 
;	common block SXR_STUFF, unless specified otherwise by the keyword ENUNS.
;    THETA
;	Angle of incidence (glancing).
; OPTIONAL INPUT PARAMETERS:
;	None.
; KEYWORD PARAMETERS:
;    ELEMENTS
;	Accepts a list of one or more elements which comprise the target.  Can
;	be provided as character array (each element represented by its 
;	chemical symbol) or numeric array (each element represented by its Z).
;    NUM_ELEMS
;	Integer scalar or array.  Contains the numbers of elements in each 
;	layer.  If only pure, single element layers are present, (i.e. if all 
;	the entries in NUM_ELEMS are 1) can be ignored.  IMPORTANT:  even if 
;	only one entry is not 1, all the entries must be provided.
;    DENS
;	Array of layer densities (same length as NUM_ELEMS).  Not needed if 
;	only pure elemental layers are present, mandatory otherwise.
;    WEIGHTS
;	Array of partial weights of the elements in each layer (same length as 
;	NUM_ELEMS), in order.  Not needed if only pure elemental layers are 
;	present, mandatory otherwise.
;    THICKS
;	Array of layer thicknesses.  Not needed if only a single interface is
;	present.
;    ENUNS
;	Character value, specifies energy units.  Acceptable units are:
;	KEV, EV, NANOMETERS, ANGSTREM.  Only first 2 letters matter.  Default
;	units are specified by the variable ENUN in the common block SXR_STUFF,
;	initially set to the value of DEFUN ('keV').  If provided, the value in
;	ENUNS replaces the previous ENUN.
;    /DEGREES
;	Switch.  Specifies that the angle is given in degrees. Default is 
;	radians.
;    /CENTIMETER
;	Switch.  Specifies that the thickness(es) are given in centimeters.  
;	Default is Angstrem.
;    /TRANSMITION
;	Switch.  Specifies a transmition mirror.  Default is reflection.
;    /COHERENT
;	Switch.  Specifies that phases add coherently between layers.  Default
;	is incoherent.  Not relevant for a single interface mirror.
; OUTPUTS:
;	Returns the reflection (or, optionally, transmition) coefficient for 
;	all the energies in EN.  Output form is same as EN.
; OPTIONAL OUTPUT PARAMETERS:
;	None.
; COMMON BLOCKS:
;       SXR_STUFF.  See LOAD_ABS_COEFFS for more information.
; SIDE EFFECTS:
;	If necessary the common block SXR_STUFF will be initialized via a call 
;       to LOAD_ABS_COEFFS.
; RESTRICTIONS:
;	The environment variable IDL_ABS_COEFFS should point to the XDR binary
;       file of absorption coefficients.  This file is read by LOAD_ABS_COEFFS.
; PROCEDURE:
;	Multiple interface backwards propagation.  
;	and absorption coefficients provided by ABS_COEFF.  Uses calls to EVUN 
;	and DIELECT.
; MODIFICATION HISTORY:
;	Created 15-MARCH-1993 by Mati Meron.
;-

    common sxr_stuff, sorlist, abctab, enun, denun, curbeam

    load_abs_coeffs
    on_error, 1
    if n_elements(abctab) eq 0 then message, 'Tables not loaded!'
    enun = Default(dun,enun,/strict)
    conv = 12.398424

    k = 2.*!pi/conv*Econv(en, from = enun, to = denun)
    if keyword_set(deg) then sang = sin(!dtor*theta) else sang = sin(theta)
    defrarr = abctab(Elecomp(elarr)).ro
    nen = n_elements(en)
    nels = n_elements(elarr) 
    nlays = n_elements(Default(numels,elarr))

    if nlays eq nels then begin
	numels = replicate(1,nlays)
	rarr = Default(rarr,replicate(0.,nlays))
	warr = Default(warr,replicate(1.,nlays))
    endif else begin
	numels = fix(numels)
	if nels ne total(numels) then message, 'Wrong # of elements!'
	if n_elements(rarr) ne nlays then message, 'Wrong # of layer densities!'
	if n_elements(warr) ne nels then message, 'Wrong # of partial weights!'
    endelse
    
    trafl = keyword_set(tran)
    nenal = nlays + trafl
    if n_elements(tarr) lt (nenal - 1) then message, 'Wrong # of thicknesses!'
    if nenal eq 1 then d = [0., 0.] else d = [0., tarr(0:nenal-2), 0.]
    if keyword_set(cm) then d = 1e8*d

    ni = make_array(nen, nenal + 1, type = 6)
    zi = ni
    ephi = ni

    jf = 0
    ni(*,0) = sang
    for i = 1, nlays do begin
	jl = jf + numels(i-1) - 1
	if jl eq jf and rarr(i-1) le 0 then begin
	    rarr(i-1) = defrarr(jf)
	    warr(jf) = 1.
	endif
	ni(*,i) = sqrt(sang^2 - Dielect(en, elements = elarr(jf:jl), $
		    dens = rarr(i-1), weights = warr(jf:jl)))
	zi(*,i) = (ni(*,i-1) - ni(*,i))/(ni(*,i-1) + ni(*,i))
	ephi(*,i) = exp(complex(0,2*d(i))*k*ni(*,i))
	jf = jl + 1
    endfor
    if trafl then begin
	ni(*,nenal) = sang
	zi(*,nenal) =(ni(*,nenal-1) - ni(*,nenal))/(ni(*,nenal-1) + ni(*,nenal))
	ephi(*,nenal) = complex(1,0)
    endif
    
    if keyword_set(coh) then begin
	r = replicate(complex(0.,0.), nen)
	t = replicate(complex(1.,0.), nen)
	for i = nenal, 1, -1 do begin
	    dumf = r*ephi(*,i)
	    r = dumf + t*zi(*,i)
	    t = dumf*zi(*,i) + t
	endfor
	r = abs(r)^2
	t = abs(t)^2
    endif else begin
	r = replicate(0., nen)
	t = replicate(1., nen)
	for i = nenal, 1, -1 do begin
	    dumf = r*abs(ephi(*,i))^2
	    dums = abs(zi(*,i))^2
	    r = dumf + t*dums
	    t = dumf*dums + t
	endfor
    endelse

    res = 1./t
    if not trafl then res = r*res else $
    for i = 1, nenal do res = res*abs(ephi(*,i)*(1. + zi(*,i))^2)

    return, res
end
