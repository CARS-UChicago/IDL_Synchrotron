Function KJK_pol, k, n, m

;+
; NAME:
;		KJK_POL
; VERSION:
;		5.1
; PURPOSE:
;		Evaluation of an expression occurring within the KJK integral.
; CATEGORY:
;		Mathematical, x-ray specific.
; CALLING SEQUENCE:
;		Result = KJK_POL( K, N [, M])
; INPUTS:
;	K
;		The K-value of an undulator
;	N
;		Power, see below.  Defaults to 0.
; OPTIONAL INPUT PARAMETERS:
;	M
;		Second optional power, see below.  Defaults to 0.
; KEYWORD PARAMETERS:
;		None.
; OUTPUTS:
;		Returns the value of the integral of (1 + K^2*cos(t)^2)^N, or, if M is
;		provided, of ((1 + K^2*cos(t)^2)^N)*cos(t)^(2*M), evaluated by residua.
; OPTIONAL OUTPUT PARAMETERS:
;		None.
; COMMON BLOCKS:
;		None.
; SIDE EFFECTS:
;		None.
; RESTRICTIONS:
;		None.
; PROCEDURE:
;		Summing a series resulting from a residua evaluation.  Calls BINCOEF,
;		CALCTYPE, CAST, DEFAULT and ISNUM from MIDL.
; MODIFICATION HISTORY:
;		Created 30-OCT-2005 by Mati Meron.
;-

	on_error, 1

	n = Default(n,0)
	typ = Calctype(0.,k)
	if n_elements(k) eq 1 then wk = k[0] else message, 'K must be a scalar!'

	if Isnum(m,/int) then begin
		if m ge 0 then begin
			if wk gt exp(alog(Toler())/(2*m)) then begin
				q = lindgen(m+1)
				tem = (-1)^(m-q)*Bincoef(1d*m,q)
				for r = 0l, m do tem[r] = tem[r]*KJK_pol(wk,1d*n+q[r])
				res = total(tem)/wk^(2*m)
			endif else res = $
			Bincoef(2d*m,m)/2.^(2*m)*(1 + wk^2*n*(2*m+1.)/(2*m+2.))
		endif else message, 'M, if defined, must be non-negative!'
	endif else begin
		zl = (1 + wk^2/2. + sqrt(1 + wk^2))/2
		zs = (wk/2.)^4/zl
		if n ge 0 then begin
			p = lindgen(n+1)
			res = total(Bincoef(1d*n,p)^2*zs^(n-p)*zl^(p))
		endif else begin
			p = lindgen(-n)
			res= total(Bincoef(1d*n,p)*Bincoef(-1d*n-1,p)*(zl-zs)^(n-p)*(-zs)^p)
		endelse
	endelse

	return, Cast(res,typ,typ,/fix)
end