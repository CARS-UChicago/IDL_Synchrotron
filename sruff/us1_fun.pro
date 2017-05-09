Function US1_fun, l, u, v, eps, verify = ver, status = st

;+
; NAME:
;		US1_FUN
; VERSION:
;		4.2
; PURPOSE:
;		Calculates the sum /Sigma{-/inf,/inf}{J_n(lu)*J_2n+l(lv) which is used
;		in SR undulator calculations.
; CATEGORY:
;		Mathematical, SR specific
; CALLING SEQUENCE:
;		Result = US1_FUN (L, U, V [, EPS] [, STATUS = ST])
; INPUTS:
;	L
;		Integer scalar.
;	U
;		Real scalar.
;	V
;		Real scalar.
; OPTIONAL INPUT PARAMETERS:
;	EPS
;		Specifies allowed relative calculation error.  Default value is set
;		according to the machine precision of U and V.
; KEYWORD PARAMETERS:
;	/VERIFY
;		Switch.  When set, passed to ROMBERG (see there) to provide for more
;		precise integration (at the cost of speed).
;	STATUS
;		Optional output, see below.
; OUTPUTS:
;		Returns calculation result, type FLOAT or the higher of [U, V].
; OPTIONAL OUTPUT PARAMETERS:
;	STATUS
;		The name of the variable to receive the calculation status.  Possible
;		values are:
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
;		Calculates the sum expressed as an integral over a finite range.  The
;		integration is performed through a call to ROMBERG (in MIDL) using the
;		kernel JJ1_ARG (in SRUFF).  Also calls CAST, DEFAULT, TOLER and TYPE
;		from MIDL.
; MODIFICATION HISTORY:
;		Created 30-MARCH-1994 by Mati Meron.
;		Modified 15-SEP-2001 by Mati Meron.  Verified WINDOWS compatibility.
;		Modified 20-JAN-2006 by Mati Meron.  Added keyword /VERIFY.
;-

	st = 1
	eps = Default(eps,(Toler(u) > Toler(v)),/dtype)
	ll = fix(l)
	if abs(v) lt eps then begin
		lh = ll/2
		if ll mod 2 eq 0 then res = (-1)^lh*beselj(ll*u,lh) $
		else res = (-1)^lh*(ll*v/2.)*(beselj(ll*u,lh) + beselj(ll*u,lh+1))
	endif else res = (l eq 0) + 1/!dpi* $
	Romberg('JJ1_arg',!dpi*[-1,1],eps,par=[ll,u,v],/rel,stat=st,try=3,ver=ver)
	return, Cast(res,4,Type([u,v]),/fix)
end