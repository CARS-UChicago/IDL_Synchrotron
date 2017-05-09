Function Pendep, theta = tet, qz = qzv, energy = ekev, wavelength = lamb, $
	elements = elar, dens = ro, dfac = dfac, weights = warr, enuns = dun, $
	degrees = deg, formula = form, ret_vals = rtv

;+
; NAME:
;		PENDEP
; VERSION:
;		5.5
; PURPOSE:
;		Calculates the penetration depth of X-rays into a flat surface.
; CATEGORY:
;		X-ray calculations.
; CALLING SEQUENCE:
;		Result = PENDEP ([THETA = TET ... QZ = QZV],
;							[ENERGY = EKEV ... WAVELENGTH = LAMB],
;							ELEMENTS = ELAR [, optional keywords])
; INPUTS:
;		None
; OPTIONAL INPUT PARAMETERS:
;		None.
; KEYWORD PARAMETERS:
;	THETA									|	Either THETA or QZ, but never
;		Angle(s) of incidence (glancing).	| 	both, must be given.  Numeric
;		In radians, unless /DEGREES is set.	|	values of arbitrary dimensions.
;	QZ										|
;		Value(s) of Qz (in inverse Angstrem)|
;
;	ENERGY									|	Either Energy or Wavelength,
;		Photon energy (in keV).				|	but never both, must be given.
;	WAVELENGTH								|	Only scalars allowed.
;		Photon wavelength (in Angstrem).	|
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
;		initially set to the value of DENUN ('keV').  If provided, the value in
;		ENUNS replaces the previous ENUN.
;	/FORMULA
;		Switch.  If set and WEIGHTS (see above) are provided, they're taken as
;		"formula weights".
;	/DEGREES
;		Switch.  Specifies that the angle is given in degrees. Default is
;		radians.
;	RET_VALS
;		Optional output, see below.
; OUTPUTS:
;		Returns the penetration depth, in Angstrem, for all the angles or qz
;		values provided.  Output form is same as input.
; OPTIONAL OUTPUT PARAMETERS:
;	RET_VALS
;		If THETA is given, returns the corresponding QZ and vice versa.
; COMMON BLOCKS:
;		None
; SIDE EFFECTS:
;		None.
; RESTRICTIONS:
;		Only elements defined in ABCTAB may be used
; PROCEDURE:
;		Straightforward, from definitions.  Calls DIELECT from SRUFF.  Calls
;		ARREQ, FPU_FIX, IMAGINARY_MM and ONE_OF from MIDL.
; MODIFICATION HISTORY:
;		Created 1-JUN-2006 by Mati Meron.
;-

	on_error, 1

	case One_of(ekev,lamb,val=input) of
		0	:	en = input
		1	:	en = conv/input
		else:	message, 'Missing energy/wavelength input!'
	endcase
	k = 2*!pi*en/!srcon.conv

	if keyword_set(deg) then amult = !dtor else amult = 1.
	case One_of(tet,qzv,val=input) of
		0	:	begin
					sang = sin(amult*input)
					rtv = 2*k*sang
				end
		1	:	begin
					sang = input/(2*k)
					rtv = asin(sang)/amult
				end
		else:	message, 'Missing angle/qz input'
	endcase

	sandim = size(sang,/dim)
	ni = sqrt(sang^2 - Dielect(en, elements = elar, dens = ro, dfac = dfac, $
		weights = warr, enuns = dun, form= form))
	res = 1/(2*k*Imaginary_mm(ni))
	if not Arreq(sandim,0) then res = reform(res,sandim)

	return, FPU_fix(res)
end