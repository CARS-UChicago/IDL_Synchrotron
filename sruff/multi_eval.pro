Function Multi_eval, theta = tet, qz = qzv, energy = ekev, wavelength = lamb, $
	degrees = deg, transmit = tran, field = fld

;+
; NAME:
;		MULTI_EVAL
; VERSION:
;		5.5
; PURPOSE:
;		Calculates the reflection of X-rays off a multilayer.
; CATEGORY:
;		X-ray calculations.
; CALLING SEQUENCE:
;		Result = MULTI_EVAL ([THETA = TET ... QZ = QZV],
;							[ENERGY = EKEV ... WAVELENGTH = LAMB],
;							[, optional keywords])
; INPUTS:
;		None
; OPTIONAL INPUT PARAMETERS:
;		None.
; KEYWORD PARAMETERS:
;	THETA									|	Either THETA or QZ, but never
;		Angle(s) of incidence (glancing).	| 	both, must be given.
;		In radians, unless /DEGREES is set.	|
;	QZ										|
;		Value(s) of Qz (in inverse Angstrem)|
;
;	ENERGY									|	Either Energy or Wavelength,
;		Photon energy (in keV).				|	but never both, must be given.
;	WAVELENGTH								|
;		Photon wavelength (in Angstrem).	|
;
;		Note:	Either THETA (or QZ) or ENERGY (or WAVELENGTH) may be vectors,
;				but not both.
;	/DEGREES
;		Switch.  Specifies that the angles are given in degrees. Default is
;		radians.
;	/TRANSMITION
;		Switch.  Specifies transmition calculation.  Default is reflection.
;	/FIELD
;		Switch.  If set, field values (either reflected or, if /TRANSMISSION is
;		set, transmitted are returned, normalized so that the incoming incident
;		field is 1.  Note, these values are complex.
; OUTPUTS:
;		Returns the reflection (or, optionally, transmition) coefficient for
;		all the values provided, unless /FIELD is set in which case field values
;		are returned.  Output form is same as larger of inputs.
; OPTIONAL OUTPUT PARAMETERS:
;		None.
; COMMON BLOCKS:
;		MULTILAYER_PARS.  See MULTI_GEN for more information.
; SIDE EFFECTS:
;		None.
; RESTRICTIONS:
;		See above about vector inputs.
; PROCEDURE:
;		Utilizes MIRROR or REFLECT from SRUFF, see there for details.  Also
;		calls CALCTYPE, CAST, FPU_FIX and from MIDL.
; MODIFICATION HISTORY:
;		Created 5-JUN-2006 by Mati Meron.
;-

	common multilayer_pars, m_exs, m_def, m_form, $
	m_elar, m_numels, m_rarr, m_warr, m_tarr, m_roff

	on_error, 1

	if keyword_set(deg) then amult = !dtor else amult = 1.

	check = One_of(ekev,lamb,val=input)
	typ = Calctype(0.,input)
	case check of
		0	:	ene = input
		1	:	ene = Cast(!srcon.conv/input,typ,typ)
		else:	message, 'Missing energy/wavelength input!'
	endcase
	nene = n_elements(ene)

	check = One_of(tet,qzv,val=input)
	typ = Calctype(0.,input)
	case check of
		0	:	ang = amult*tet
		1	:	if nene eq 1 then $
				ang = asin(Cast(!srcon.conv*input/(4*!dpi*ene),typ,typ)) $
				else message, 'Unacceptable option!'
		else:	message, 'Missing angle/qz input'
	endcase
	nang = n_elements(ang)

	if nene eq 1 or nang eq 1 then begin
		if nang gt 1 then res = Reflect(theta=ang,ene=ene, $
		elem=m_elar,num=m_numels,dens=m_rarr,wei=m_warr,thi=m_tarr,rou=m_roff,$
		form=m_form,tran=tran,field=fld) $
		else res = Mirror(ene,ang, $
		elem=m_elar,num=m_numels,dens=m_rarr,wei=m_warr,thi=m_tarr,rou=m_roff,$
		form=m_form,tran=tran,field=fld)
	endif else message, 'Either energy or angle must be scalars!'

	return, FPU_fix(res)
end