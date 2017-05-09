Function H_fun, gpsi, e_ec, pow, normalized = norm

;+
; NAME:
;		H_FUN
; VERSION:
;		6.0
; PURPOSE:
;		Calculates the bending magnet power distribution function H (scaled by
;		(E/(2*E_c))^2), from R.Dejus at. al. paper, NIM A319 (1992) 207-212.
; CATEGORY:
;		Mathematical, SR specific
; CALLING SEQUENCE:
;		Result = H_FUN (GPSI, E_EC)
; INPUTS:
;	GPSI													|
;		Scalar or vector.  /Gamma*/Psi						|	At least one of
;	E_EC													|	these two must
;		Scalar or vector.  Photon energy in units of the 	|	be a scalar.
;		critical energy E_c.								|
;	POW
;		Scalar.  If provided, the result is raised to power POW prior to
;		returning.  For internal use.
; OPTIONAL INPUT PARAMETERS:
;		None.
; KEYWORD PARAMETERS:
;	/NORMALIZED
;		Switch.  If set, the result is normalized to yield 1 for GPSI = 0.
;		Only active if E_EC is a scalar.
; OUTPUTS:
;		Returns calculation result, type FLOAT or the highest of the types of
;		[GPSI, E_EC].
; OPTIONAL OUTPUT PARAMETERS:
;		None
; COMMON BLOCKS:
;		None.
; SIDE EFFECTS:
;		None.
; RESTRICTIONS:
;		None.
; PROCEDURE:
;		Straighforward, follows the formula.  Uses ARREQ, BESELK_MM, CAST,
;		FPU_FIX and ISNUM from MIDL.
; MODIFICATION HISTORY:
;		Created 25-JAN-1994 by Mati Meron.
;		Modified 15-SEP-2001 by Mati Meron.  Verified WINDOWS compatibility.
;		Modified 10-JAN-2005 by Mati Meron.  Reversed order of inputs for the
;		purpose of integration.  Added variable POW.
;		Modified 20-JAN-2007 by Mati Meron.  Internal changes.
;-

	on_error, 1

	asc = Arreq((size(gpsi))[[0,2]],[0,1])
	esc = Arreq((size(e_ec))[[0,2]],[0,1])
	if not (asc or esc) then message, 'Either Energy or Gam*psi must be scalar!'

	tem = Cast(1. + gpsi^2, 4)
	ksi_0 = 0.5*e_ec
	ksi = ksi_0*tem^(3./2)

	if asc then begin
		res = ksi_0^2
		dum = where(res gt 0, ndum)
		if ndum gt 0 then res[dum] = res[dum]*tem*$
		(tem*Beselk_mm(ksi[dum],2./3.)^2 + (tem-1)*Beselk_mm(ksi[dum],1./3.)^2)
	endif else begin
		if ksi_0 gt 0 then begin
			res = tem*$
			(tem*Beselk_mm(ksi,2./3.)^2 + (tem-1)*Beselk_mm(ksi,1./3.)^2)
			if keyword_set(norm) then res = res/Beselk_mm(ksi_0,2./3.)^2 $
			else res = ksi_0^2*res
		endif else res = 0*tem
	endelse
	if Isnum(pow) then if pow ne 1 then res = res^pow

	return, FPU_fix(res)
end