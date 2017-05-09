Function Diff_scatt, en, theta, rough, enuns = dun, $
    degrees = deg, nanometer = nm

;+
; NAME:
;		DIFF_SCATT
; VERSION:
;		4.2
; PURPOSE:
;		Calculates the diffuse scattering correction to reflectivity.
; CATEGORY:
;		X-ray calculations.
; CALLING SEQUENCE:
;		Result = DIFF_SCATT (EN, THETA, ROUGH [, optional keywords])
; INPUTS:
;	EN
;		Energy, assumed in the units specified by the variable ENUN in the
;		common block SXR_STUFF, unless specified otherwise by the keyword ENUNS.
;	THETA
;		Angle of incidence (glancing).
;	ROUGH
;		RMS surface roughness, assumed in Angstrem unless the keyword
;		NANOMETER is set.
; KEYWORD PARAMETERS:
;	ENUNS
;		Character value, specifies energy units.  Acceptable units are:
;		KEV, EV, NANOMETERS, ANGSTREM.  Only first 2 letters matter.  Default
;		units are specified by the variable ENUN in the common block SXR_STUFF,
;		initially set to the value of DEFUN ('keV').  If provided, the value in
;		ENUNS replaces the previous ENUN.
;	/DEGREES
;		Switch.  If set degree units are assumed for angle. Default is radians.
;	/NANOMETER
;		Switch.  Specifies that the roughness is given in nanometers.
;		Default is Angstrem.
; OUTPUTS:
;		Returns the roughness multiplier for the reflection coefficient, for
;		all the energies in EN.  Output form is same as EN.
; OPTIONAL OUTPUT PARAMETERS:
;		None.
; COMMON BLOCKS:
;		SXR_STUFF.  See LOAD_ABS_COEFFS for more information.
; SIDE EFFECTS:
;		None.
; RESTRICTIONS:
;		The common block SXR_STUFF must be initialized through a call to
;		LOAD_ABS_COEFFS before MIRROR may be used.
; PROCEDURE:
;		Straightforward from definition.  Calls ECONV and LOAD_ABS_COEFFS.
;		Also calls DEFAULT from MIDL.
; MODIFICATION HISTORY:
;		Created 10-JAN-1997 by Mati Meron.
;		Modified 15-SEP-2001 by Mati Meron.  Automated the call to
;		LOAD_ABS_COEFFS and verified WINDOWS compatibility.
;-

	common sxr_stuff, sorlist, abctab, enun, denun, curbeam

	on_error, 1
	Load_abs_coeffs
	enun = Default(dun,enun,/strict)
	conv = 12.398424

	k = 2.*!pi/conv*Econv(en, from = enun, to = denun)
	if keyword_set(deg) then sang = sin(!dtor*theta) else sang = sin(theta)
	if keyword_set(nm) then wrough = 0.1*rough else wrough = rough

	return, exp(-2*wrough*sang*k)

end