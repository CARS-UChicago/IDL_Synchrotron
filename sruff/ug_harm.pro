Function UG_harm, l, u, v, k, gtet, eps, verify = ver, status = st

;+
; NAME:
;		UG_HARM
; VERSION:
;		4.2
; PURPOSE:
;		Calculates the function G_k (see "Handbook of synchrotron radiation",
;		Vol. 1A, p. 150) which is used in SR undulator calculations.
; CATEGORY:
;		Mathematical, SR specific
; CALLING SEQUENCE:
;		Result = UG_HARM (L, U, V, K, GTET [, EPS] [, STATUS = ST])
; INPUTS:
;	L
;		Integer scalar, harmonic number.
;	U
;		Real scalar, corresponds to /ksi_z.
;	V
;		Real scalar, corresponds to /ksi_x.
;	K
;		Real scalar, undulator K.
;	GTET
;		Real scalar, equals /gamma*/theta.
; OPTIONAL INPUT PARAMETERS:
;	EPS
;		Specifies allowed relative calculation error.  Default value is
;		established according to the machine precision of U and V.
;	/VERIFY
;		Switch.  When set, passed to US1_FUN and US2_FUN to provide for more
;		precise integration (at the cost of speed).
; KEYWORD PARAMETERS:
;	STATUS
;		Optional output, see below.
; OUTPUTS:
;		Returns calculation result, type FLOAT or the higher of [U, V].
; OPTIONAL OUTPUT PARAMETERS:
;	STATUS
;		The name of the variable to receive the calculation status.  Returns
;		a two-element vector where the first value is the status of US1 and the
;		second value is the status of US2.  Possible values are:
;			0 -	Calculation didn't converge.
;			1 - OK.
;			2 - Calculation converged but with precision worse than specified.
; COMMON BLOCKS:
;		None.
; SIDE EFFECTS:
;		None.
; RESTRICTIONS:
;		None.
; PROCEDURE:
;		Calculates the G_k function using calls to US1_FUN and US2_FUN (both
;		in SRUFF).  Also calls DEFAULT, FPU_FIX and TOLER from MIDL.
; MODIFICATION HISTORY:
;		Created 30-MARCH-1994 by Mati Meron.
;		Modified 15-SEP-2001 by Mati Meron.  Verified WINDOWS compatibility.
;		Modified 20-JAN-2006 by Mati Meron.  Added keyword /VERIFY.
;-

	on_error, 1
	st1 = 1
	st2 = 1
	eps = Default(eps,(Toler(u) > Toler(v)),/dtype)
	if k eq 0 then coef = 0 else coef = 8.*l*u/k
	if coef ne 0 then begin
		fir = US1_fun(l,u,v,eps,ver=ver,stat=st1)
		ll = fix(l)
		if abs(v) lt eps then begin
			lh = ll/2
			if ll mod 2 eq 0 then sec = 0. else $
			sec = (-1)^(lh)/(2.)*(beselj(ll*u,lh) - beselj(ll*u,lh+1))
		endif else sec = (fir + 2./ll*US2_fun(l,u,v,eps,ver=ver,stat=st2))/v
		res = coef^2*((fir*gtet/k)^2 - fir*sec*v/(4.*u) + sec^2) > 0.
	endif else res = 0.
	st = [st1,st2]

	return, FPU_fix(res)
end