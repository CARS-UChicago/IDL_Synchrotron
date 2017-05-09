Function Disper_fun, x, pow, warn = wrn

;+
; NAME:
;		DISPER_FUN
; VERSION:
;		4.2
; PURPOSE:
;		Calculates the dispertion integral part of the dielectric coefficient.
; CATEGORY:
;		Mathematical, SR specific
; CALLING SEQUENCE:
;		Result = DISPER_FUN( X, POW)
; INPUTS:
;	X
;		Numeric, otherwise arbitrary.
;	POW
;		Scalar.  The power exponent that's used in the calculation which
;		assumes an absorption coefficient proportional to Energy^(-2 - POW).
; OPTIONAL INPUT PARAMETERS:
;		None
; KEYWORD PARAMETERS:
;	/WARN
;		Switch.  If set, a warning is issued when POW is out of the legitimate
;		(-1,1) range.
; OUTPUTS:
;		Returns an approximation to the integral from X to infinity of
;		1/(t^pow*(1-t^2)), one value for each entry in X.  If POW is out of
;		range, returns 0.
; OPTIONAL OUTPUT PARAMETERS:
;		None
; COMMON BLOCKS:
;		None.
; SIDE EFFECTS:
;		None.
; RESTRICTIONS:
;		POW must be in the range (-1, 1).  X should be positive.
; PROCEDURE:
;		Utilizing an analytical approximation to the integral above.  Calls
;		itself recursively (sometimes).  Uses CAST, TOLER and TYPE from MIDL.
; MODIFICATION HISTORY:
;		Created 30-MARCH-1993 by Mati Meron.
;		Rewritten and upgraded, 15-JUNE-1994 by Mati Meron.
;		Modified 15-SEP-2001 by Mati Meron.  Verified WINDOWS compatibility.
;		Modified 25-JUL-2002 by Mati Meron.  Added keyword WARN.
;-

	on_error, 1

	if min(x) lt 0 then message, 'X must be nonnegative!'
	sml = Toler(1d)

	if abs(pow) lt 1 then begin
		if pow ge 0 then begin
			dx = double(x)
			res = 0.5*alog(abs((1-dx)/(1+dx)) > sml)
			if pow gt 0 then begin
				fac = (1 - 2*pow)/pow
				if fac eq 0 then res = res + 1/(1 + dx) else $
				res = res + alog(abs(1 + fac/(1 + dx)) > sml)/fac
			endif
		endif else res = !pi/2*tan(!pi/2*pow) + Disper_fun(1./(x > sml), -pow)
	endif else begin
		res = 0.*x
		if keyword_set(wrn) then message, 'Parameter out of range!', /continue
	endelse

	return, Cast(res,4,Type(x),/fix)
end