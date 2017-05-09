Function Wig_fun, gtetok, gpsi, e_ec, pow, normalized = norm

;+
; NAME:
;		WIG_FUN
; VERSION:
;		6.0
; PURPOSE:
;		Calculates the wiggler power distribution function H (scaled by
;		(E/(2*E_c))^2), from R.Dejus at. al. paper, NIM A319 (1992) 207-212.
; CATEGORY:
;		Mathematical, SR specific
; CALLING SEQUENCE:
;		Result = WIG_FUN (GTETOK, GPSI, E_EC [. keywords])
; INPUTS:
;	GTETOT												|	Either E_EC is a
;		Scalar or array.  Gamma*Theta/K					|	scalar and the other
;	GPSI												|	two are arrays of
;		Scalar or vector.  Gamma*Psi					|	same structure, or
;	E_EC												|	E_EC is a vector and
;		Scalar or vector.  Photon energy in units of 	|	the other two are
;		the	critical energy E_c.						|	scalars.
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
;		[GTETOK, GPSI, E_EC].
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
;		Created 25-JAN-2007 by Mati Meron.
;-

	on_error, 1

	if (size(gtetok))[0] eq 0 then wtet = gtetok else wtet = reform(gtetok)
	if (size(gpsi))[0] eq 0 then wpsi = gpsi else wpsi = reform(gpsi)
	if Arreq((size(e_ec))[[0,2]],[0,1]) then begin
		if Arreq(wtet,wpsi,/noval) then asc = 0 $
		else message, 'Psi and Theta must have same structure!'
	endif else begin
		tetsc = Arreq((size(wtet))[[0,2]],[0,1])
		psisc = Arreq((size(wpsi))[[0,2]],[0,1])
		if (tetsc and psisc) then asc = 1 $
		else message, 'For vector energy, angles must be scalars!'
	endelse

	tem = Cast(1. + wpsi^2, 4)
	if asc then begin
		if abs(wtet) lt 1 then begin
			ksi_0 = 0.5*e_ec/sqrt(1 - wtet^2)
			ksi = ksi_0*tem^(3./2)
			res = ksi_0^2
			dum = where(res gt 0, ndum)
			if ndum gt 0 then res[dum] = res[dum]*tem*$
			(tem*Beselk_mm(ksi[dum],2./3.)^2 + $
			(tem-1)*Beselk_mm(ksi[dum],1./3.)^2)
		endif else res = 0.*e_ec
	endif else begin
		if e_ec gt 0 then begin
			dum = where(abs(wtet) lt 1, ndum)
			if ndum gt 0 then begin
				res = 0.*tem
				ksi_0 = 0.5*e_ec/sqrt(1 - wtet[dum]^2)
				ksi = ksi_0*tem[dum]^(3./2)
				res[dum] = tem[dum]*$
				(tem[dum]*Beselk_mm(ksi,2./3.)^2 + $
				(tem[dum]-1)*Beselk_mm(ksi,1./3.)^2)
				if keyword_set(norm) then $
				res[dum] = res[dum]/Beselk_mm(ksi_0,2./3.)^2 $
				else res[dum] = ksi_0^2*res[dum]
			endif
		endif else res = 0.*tem
	endelse
	if Isnum(pow) then if pow ne 1 then res = res^pow

	return, FPU_fix(res)
end