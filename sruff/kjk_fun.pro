Function KJK_fun, k, gtet, gpsi, eps, status = st, correct = cor

;+
; NAME:
;	KJK_FUN
; PURPOSE:
;	Calculates the undulator power distribution function f_k, from 
;	K.J.Kim's paper NIM A246 (1986) 67-70.
; CATEGORY:
;	Mathematical, SR specific
; CALLING SEQUENCE:
;	Result = KJK_FUN (K, GTET, GPSI, [, EPS] [ keywords])
; INPUTS:
;    K
;	Undulator K.
;    GTET
;	/Gamma*/Theta
;    GPSI
;	/Gamma*/Psi
; OPTIONAL INPUT PARAMETERS:
;    EPS
;	Specifies allowed relative calculation error.  Default value is ~2e-8.
; KEYWORD PARAMETERS:
;    STATUS
;	Optional output, see below.
;    /CORRECT
;	Specifies whether the quasi-singular part of the integral should be 
;	calculated separately.  Default is NO CORRECT.  Ignored if K<GTET.
; OUTPUTS:
;	Returns calculation result, type FLOAT or the highest of the types of
;	[K,GTET,GPSI].
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
;	Performs the integration through call to ROMBERG (in MIDL) using the
;	kernel KJK_ARG (in SRUFF).  If /CORRECT is set, calculates separately
;	the quasi-singular part using NATAN (from MIDL).  Also calls CAST,
;	DEFAULT and TYPE from MIDL.
; MODIFICATION HISTORY:
;	Created 30-MARCH-1994 by Mati Meron.
;-

    st = 1
    deps = 0.5^(51./2.)
    eps = Default(eps,deps,/dtype)
    pars = abs(Cast([k,gtet,gpsi],4))
    if keyword_set(cor) and pars(0) gt max([0,pars(1)]) then $
	pars = [pars,1] else pars = [pars,0]
    ks = 1. + pars(0)^2
    coef = 32*ks^(7./2.)/(!dpi*(5+ks*(1+ks*(3+ks*7))))
    res = Romberg('kjk_arg',!dpi*[0,1],eps/coef,par=pars,stat=st,try=3)

    if pars(3) then begin
	tfac = pars(0)^2 - pars(1)^2
	pfac = 1./(1 + pars(2)^2)
	frat = sqrt(tfac*pfac)
	clims = sqrt(pfac)*(pars(0)*[-1,1] + pars(1))
	c2 = Natan(clims,2)
	c34 = Natan(clims,3) - Natan(clims,4)
	corr = (pfac/pars(0))^2*frat*(c2(1) - c2(0) - 4*pfac*(c34(1) - c34(0)))
	res = res + corr
    endif

    return, Cast(coef*res,4,Type(pars))
end
