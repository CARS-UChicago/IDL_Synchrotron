Function Disper_fun, x, pow

;+
; NAME:
;	DISPER_FUN
; PURPOSE:
;	Calculates the dispertion integral part of the dielectric coefficient.
; CATEGORY:
;	Mathematical, SR specific
; CALLING SEQUENCE:
;	Result = DISPER_FUN( X, POW)
; INPUTS:
;    X
;	Numeric, otherwise arbitrary.
;    POW
;	Scalar.  The power exponent that's used in the calculation which 
;	assumes an absorption coefficient proportional to Energy^(-2 - POW).
; OPTIONAL INPUT PARAMETERS:
;	None
; KEYWORD PARAMETERS:
;	None.
; OUTPUTS:
;	Returns an approximation to the integral from X to infinity of
;	1/(t^pow*(1-t^2)), one value for each entry in X.
; OPTIONAL OUTPUT PARAMETERS:
;	None
; COMMON BLOCKS:
;	None.
; SIDE EFFECTS:
;	None.
; RESTRICTIONS:
;	POW must be in the range (-1, 1).  X should be positive.
; PROCEDURE:
;	Utilizing an analytical approximation to the integral above.  Calls 
;	itself recursively (sometimes).  Uses CAST and TYPE from MIDL.
; MODIFICATION HISTORY:
;	Created 30-MARCH-1993 by Mati Meron.
;	Rewritten and upgraded, 15-JUNE-1994 by Mati Meron.
;-

    if abs(pow) ge 1 then message, 'Parameter out of range!'
    if min(x) lt 0 then message, 'X must be nonnegative!'
    sml = 0.5^54

    if pow ge 0 then begin
	dx = double(x)
	res = 0.5*alog(abs((1-dx)/(1+dx)) > sml)
	if pow gt 0 then begin
	    fac = (1 - 2*pow)/pow
	    if fac eq 0 then res = res + 1/(1 + dx) else $
	    res = res + alog(abs(1 + fac/(1 + dx)) > sml)/fac
	endif
    endif else res = !pi/2*tan(!pi/2*pow) + Disper_fun(1./(x > sml), -pow)

    return, Cast(res,4,Type(x))
end
