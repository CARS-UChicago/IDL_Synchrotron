Function US1_fun, l, u, v, eps, status = st

;+
; NAME:
;	US1_FUN
; PURPOSE:
;	Calculates the sum /Sigma{-/inf,/inf}{J_n(lu)*J_2n+l(lv) which is used
;	in SR undulator calculations.
; CATEGORY:
;	Mathematical, SR specific
; CALLING SEQUENCE:
;	Result = US1_FUN (L, U, V [, EPS] [, STATUS = ST])
; INPUTS:
;    L
;	Integer scalar.
;    U
;	Real scalar.
;    V
;	Real scalar.
; OPTIONAL INPUT PARAMETERS:
;    EPS
;	Specifies allowed relative calculation error.  Default value is ~2e-8.
; KEYWORD PARAMETERS:
;    STATUS
;	Optional output, see below.
; OUTPUTS:
;	Returns calculation result, type FLOAT or the higher of [U, V].
; OPTIONAL OUTPUT PARAMETERS:
;    STATUS
;	The name of the variable to receive the calculation status.  Possible
;	values are:
;	    0 -	Calculation didn't converge.
;	    1 - OK.
;	    2 - Calculation converged but with precision worse than specified.
; COMMON BLOCKS:
;	None.
; SIDE EFFECTS:
;	None.
; RESTRICTIONS:
;	None.
; PROCEDURE:
;	Calculates the sum expressed as an integral over a finite range.  The
;	integration is performed through a call to ROMBERG (in MIDL) using the
;	kernel JJ1_ARG (in SRUFF).  Also calls CAST, DEFAULT and TYPE from MIDL.
; MODIFICATION HISTORY:
;	Created 30-MARCH-1994 by Mati Meron.
;-

    st = 1
    deps = 0.5^(51./2.)
    eps = Default(eps,deps,/dtype)
    if abs(v) lt deps then begin
	lh = fix(l)/2
	if fix(l) mod 2 eq 0 then res = (-1)^lh*beselj(l*u,lh) else res = 0.
    endif else res = (l eq 0) + $
    1/!dpi*Romberg('JJ1_arg',!dpi*[-1,1],eps,par=[l,u,v],/rel,stat=st,try=3)
    return, Cast(res,4,Type([u,v]))
end
