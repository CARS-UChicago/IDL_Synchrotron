Function H_fun, e_ec, gpsi, normalized = norm

;+
; NAME:
;	H_FUN
; PURPOSE:
;	Calculates the bending magnet power distribution function H, from 
;	R.Dejus at. al. paper, NIM A319 (1992) 207-212.
; CATEGORY:
;	Mathematical, SR specific
; CALLING SEQUENCE:
;	Result = H_FUN (E_EC, GPSI)
; INPUTS:
;    E_EC
;	Scalar.  Photon energy in units of the critical energy E_c.
;    GPSI
;	Scalar or vector.  /Gamma*/Psi
; OPTIONAL INPUT PARAMETERS:
;	None.
; KEYWORD PARAMETERS:
;    /NORMALIZED
;	Switch.  If set, the result is normalized to yield 1 for GPSI = 0.
; OUTPUTS:
;	Returns calculation result, type FLOAT or the highest of the types of
;	[E_EC, GPSI].
; OPTIONAL OUTPUT PARAMETERS:
;	None
; COMMON BLOCKS:
;	None.
; SIDE EFFECTS:
;	None.
; RESTRICTIONS:
;	None.
; PROCEDURE:
;	Straighforward, follows the formula.  Uses CAST from MIDL and BESELK 
;	from SRUFF.
; MODIFICATION HISTORY:
;	Created 25-JAN-1994 by Mati Meron.
;-

    on_error, 1
    if n_elements(e_ec) ne 1 then message, 'Energy must be scalar!'

    tem = Cast(1. + gpsi^2, 4)
    ksi = 0.5*e_ec*tem^(3./2)
    if keyword_set(norm) then nof = 1./Beselk(0.5*e_ec,2./3.)^2 else nof = 1.

    return, nof*tem*(tem*Beselk(ksi,2./3.)^2 + (tem - 1)*Beselk(ksi,1./3.)^2)
end
