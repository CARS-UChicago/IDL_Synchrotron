Function Mirror, en, theta, enuns = dun, elements = elar, num_elems = numels, $
    dens = rarr, dfac = dfac, weights = garr, thicks = tarr, roughness = roff, $
    multi = mlt, formula = form, degrees= deg, centimeter= cm, transmit= tran, $
    field = fld
;+
; NAME:
;		MIRROR
; VERSION:
;		5.5
; PURPOSE:
;		Calculates the reflection (or transmition) coefficient of an X-ray
;		mirror.
; CATEGORY:
;		X-ray calculations.
; CALLING SEQUENCE:
;		Result = MIRROR (EN, THETA, ELEMENTS = ELAR [, optional keywords])
; INPUTS:
;	EN
;		Energy, assumed in the units specified by the variable ENUN in the
;		common block SXR_STUFF, unless specified otherwise by the keyword ENUNS.
;	THETA
;		Angle of incidence (glancing).
; OPTIONAL INPUT PARAMETERS:
;		None.
; KEYWORD PARAMETERS:
;	ENUNS
;		Character value, specifies energy units.  Acceptable units are:
;		KEV, EV, NANOMETERS, ANGSTREM.  Only first 2 letters matter.  Default
;		units are specified by the variable ENUN in the common block SXR_STUFF,
;		initially set to the value of DENUN ('keV').  If provided, the value in
;		ENUNS replaces the previous ENUN.
;	ELEMENTS
;		Accepts a list of one or more elements which comprise the target.  Can
;		be provided as character array (each element represented by its
;		chemical symbol) or numeric array (each element represented by its Z).
;	NUM_ELEMS
;		Integer scalar or array.  Contains the numbers of elements in each
;		layer.  If only pure, single element layers are present, (i.e. if all
;		the entries in NUM_ELEMS are 1) can be ignored.  IMPORTANT:  even if
;		only one entry is not 1, all the entries must be provided.
;	DENS
;		Array of layer densities (same length as NUM_ELEMS).  Not needed if
;		only pure elemental layers are present, mandatory otherwise.
;	DFAC
;		Density multiplier, used to account for the possibility that layer
;		densities differ from bulk densities.  Can be either a scalar or an
;		array of same length as NUM_ELEMS.  Default value is 1.
;	WEIGHTS
;		Array of partial weights of the elements in each layer (same length as
;		NUM_ELEMS), in order.  Not needed if only pure elemental layers are
;		present, mandatory otherwise.
;	THICKS
;		Array of layer thicknesses, in Angstrem (unless /CENTIMETER is set).
;		Not needed if only a single interface is present.
;	ROUGHNESS
;		Value of mirror roughness.  Assumed in Angstrem unless NANOMETERS is
;		used for ENUNS (see below).  Default value is 0.  Can be provided as
;		scalar (in which case the value provided applies to all interfaces or
;		a vector of same length as number of interfaces (in which case
;		consecutive roughnesses are applied to consecutive interfaces).
;	MULT
;		If provided and greater than 1, the target is iterated MULT times, thus
;		providing a multilayer result.  Default value is 1.
;	/FORMULA
;		Switch.  If set and WEIGHTS (see above) are provided, they're taken as
;		"formula weights".
;	/DEGREES
;		Switch.  Specifies that the angle is given in degrees. Default is
;		radians.
;	/CENTIMETER
;		Switch.  Specifies that the thickness(es) are given in centimeters.
;		Default is Angstrem.
;	/TRANSMIT
;		Switch.  Specifies a transmition mirror.  Default is reflection.
;	/FIELD
;		Switch.  If set, field values (either reflected or, if /TRANSMISSION is
;		set, transmitted are returned, normalized so that the incoming incident
;		field is 1.  Note, these values are complex.
; OUTPUTS:
;		Returns the reflection (or, optionally, transmition) coefficient for
;		all the energies in EN, unless /FIELD is set in which case field values
;		are returned.  Output form is same as EN.
; OPTIONAL OUTPUT PARAMETERS:
;		None.
; COMMON BLOCKS:
;		SXR_STUFF.  See LOAD_ABS_COEFFS for more information.
; SIDE EFFECTS:
;		None.
; RESTRICTIONS:
;		Only elements defined in ABCTAB may be used.
; PROCEDURE:
;		Multiple interface backwards propagation.  Uses values from ABCTAB and
;		absorption coefficients provided by ABS_COEFF.  Uses calls to ECONV,
;		ELECOMP, DIELECT and LOAD_ABS_COEFFS.  Also calls ABS_MM, ARREQ, CAST,
;		DEFAULT, FPU_FIX, ISNUM, REAL_MM and STREQ from MIDL.
; MODIFICATION HISTORY:
;		Created 15-MAR-1993 by Mati Meron.
;		Modified 15-JUN-1996 by Mati Meron.  Changed default from COHERENT to
;		INCOHERENT.
;		Modified 15-SEP-2001 by Mati Meron.  Verified WINDOWS compatibility.
;		Modified 20-JUL-2002 by Mati Meron.  Added roughness on *all*
;		interfaces capability and eliminated the keyword INCOHERENT which is
;		no longer needed.
;		Modified 25-JUL-2002 by Mati Meron.  Added keyword FORMULA.
;		Modified 5-JUN-2003 by Mati Meron.  Added keyword DFAC.
;		Modified 15-JUN-2003 by Mati Meron.  Added keyword MULT.
;		Modified 1-JUN-2006 by Mati Meron.  Added keyword FIELD.
;-

	common sxr_stuff, sorlist, abctab, enun, denun, curbeam

	on_error, 1
	Load_abs_coeffs
	enun = Default(dun,enun,/strict)
	conv = 12.398424
	in_units = enun
	out_units=denun

	k = 2.*!pi/conv*Econv(en, from = in_units, to = out_units)
	if keyword_set(deg) then sang = sin(!dtor*theta) else sang = sin(theta)
	roffl = Isnum(roff)
	trafl = keyword_set(tran)
	mlt = Default(mlt,1l,/dtyp) > 1
	mltfl = mlt gt 1

	endim = size(en,/dim)
	nen = n_elements(en)
	nels = n_elements(elar)
	nlays = n_elements(Default(numels,elar))
	defrarr = abctab[Elecomp(elar)].ro

	if nlays eq nels then begin
		numels = Default(numels,replicate(1l,nlays),/dtyp)
		wrar = Default(rarr,replicate(0.,nlays),/dtyp)
		wgar = Default(garr,replicate(1.,nels),/dtyp)
	endif else begin
		numels = long(numels)
		wrar = Cast(rarr,4)
		wgar = Cast(garr,4)
	endelse
	if total(numels) ne nels then message, 'Wrong # of elements!'
	if n_elements(wrar) ne nlays then message, 'Wrong # of layer densities!'
	if n_elements(wgar) ne nels then message, 'Wrong # of partial weights!'

	case n_elements(dfac) of
		0	:	wdfc = replicate(1.,nlays)
		1	:	wdfc = replicate(dfac,nlays)
		nlays:	wdfc = dfac
		else:	message, 'Wrong # of density factors!'
	endcase

	interf = nlays + (trafl or mltfl)
	if n_elements(tarr) lt interf - 1 then message, 'Wrong # of thicknesses!'
	if n_elements(tarr) eq 0 then d = [0.,0.] else d = ([0.,tarr,0.])[0:nlays]
	if keyword_set(cm) then d = 1e8*d

	nenal = nlays + trafl
	ni = (zi = (ephi = (etphi = make_array(nen, nenal + 1, type = 6))))
	if roffl then begin
		case n_elements(roff) of
			1			:	wrof = [0.,replicate(roff,interf)]
			interf - 1	:	wrof = [0.,roff,roff[0]]
			interf		:	wrof = [0.,roff]
			else		:	message, 'Wrong # of roughnesses'
		endcase
		if Streq(enun,'Nano',2) then wrof = 10*wrof
		chi = (psi = ni)
	endif else chi = (psi = make_array(nen,nenal + 1,val = 1.))

	jf = 0l
	ni[*,0] = sang
	for i = 1l, nlays do begin
		jl = jf + numels[i-1] - 1
		if jl eq jf and wrar[i-1] le 0 then begin
			wrar[i-1] = defrarr[jf]
			wgar[jf] = 1.
		endif
		ni[*,i] = sqrt(sang^2 - Dielect(en, elements = elar[jf:jl], $
		dens = wrar[i-1], dfac = wdfc[i-1], weights = wgar[jf:jl], form= form))
		zi[*,i] = (ni[*,i-1] - ni[*,i])/(ni[*,i-1] + ni[*,i])
		ephi[*,i] = exp(complex(0,k*d[i])*ni[*,i])
		etphi[*,i] = (ephi[*,i])^2
		if roffl then begin
			chi[*,i] = exp(-(2*k*wrof[i]*Real_mm(ni[*,i-1]))^2/2)
			psi[*,i] = exp(-(k*wrof[i]*Real_mm(ni[*,i] - ni[*,i-1]))^2/2)
		endif
		jf = jl + 1
	endfor
	if mltfl then begin
		zi[*,0] = zi[*,1]
		zi[*,1] = (ni[*,nlays] - ni[*,1])/(ni[*,nlays] + ni[*,1])
		if roffl then begin
			chi[*,0] = chi[*,1]
			psi[*,0] = psi[*,1]
			chi[*,1] = exp(-(2*k*wrof[nlays+1]*Real_mm(ni[*,nlays]))^2/2)
			psi[*,1]= exp(-(k*wrof[nlays+1]*Real_mm(ni[*,1] - ni[*,nlays]))^2/2)
		endif
		ephi[*,0] = ephi[*,1]
		etphi[*,0] = etphi[*,1]
	endif
	rot = replicate(complex(0.,0.), nen)
	if trafl then begin
		i = nenal
		ni[*,i] = sang
		zi[*,i] =(ni[*,i-1] - ni[*,i])/(ni[*,i-1] + ni[*,i])
		ephi[*,i] = (etphi[*,i] = complex(1,0))
		if roffl then begin
			chi[*,i] = exp(-(2*k*wrof[i]*Real_mm(ni[*,i-1]))^2/2)
			psi[*,i] = exp(-(k*wrof[i]*Real_mm(ni[*,i] - ni[*,i-1]))^2/2)
		endif
		oot = replicate(complex(1.,0.), nen)
	endif

	nm = nlays*mlt
	tem = lindgen(nlays) + 1
	j = lonarr(nm + trafl)
	for i = 0l, mlt - 1 do j[i*nlays:(i+1)*nlays-1] = tem
	if trafl then j[nm] = nlays + 1
	if mltfl then j[0] = 0
	for i = nm + trafl - 1, 0, -1 do begin
		dumf = psi[*,j[i]]*etphi[*,j[i]]*rot
		dums = zi[*,j[i]]*dumf + 1
		rot = chi[*,j[i]]*(dumf + zi[*,j[i]])/dums
		if trafl then oot = psi[*,j[i]]*(zi[*,j[i]] + 1)*ephi[*,j[i]]/dums*oot
	endfor

	if trafl then res = oot else res = rot
	if not keyword_set(fld) then res = Abs_mm(res)^2
	if not Arreq(endim,0) then res = reform(res,endim)

	return, FPU_fix(res)
end