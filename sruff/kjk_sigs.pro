Function KJK_sigs, k, corrected = cor

;+
; NAME:
;		KJK_SIGS
; VERSION:
;		5.1
; PURPOSE:
;		Evaluation of angular "sigma" values for an undulator.
; CATEGORY:
;		Mathematical, x-ray specific.
; CALLING SEQUENCE:
;		Result = KJK_SIGS( K )
; INPUTS:
;	K
; OPTIONAL INPUT PARAMETERS:
;		None.
; KEYWORD PARAMETERS:
;	/CORRECTED
;		Switch.  If set, the calculated sigma values are adjusted so as to yield
;		nearly correct integrals.
; OUTPUTS:
;		Returns a two-element vector including, in order:
;			1) Gamma*Sigma(theta_x) (horizontal)
;			2) Gamma*Sigma(theta_y) (vertical)
;		for the ID radiation pattern (represented as a gaussian).
; OPTIONAL OUTPUT PARAMETERS:
;		None.
; COMMON BLOCKS:
;		None.
; SIDE EFFECTS:
;		None.
; RESTRICTIONS:
;		None.
; PROCEDURE:
;		Calculates the Sigma values based on the second derivatives of the
;		undulator power distribution function f_k, from K.J.Kim's paper
;		NIM A246 (1986) 67-70.  If /CORRECTED is set, scales the calculated
;		values so as to bring the calculated integrals (over a finite aperture)
;		into closer agreement with exact calculations.  Calls KJK_POL.  Also
;		calls CALCTYPE, and VINP from MIDL.
; MODIFICATION HISTORY:
;		Created 10-NOV-2005 by Mati Meron.
;-

	on_error, 1

	cnr = [1,-4,4,0,0]
	sdx = [0,42,-336,760,-480]
	sdy = [0,-6,40,-40,0]

	hfn = dblarr(5)
	for i = 0, 4 do hfn[i] = KJK_pol(1d*k,-(i+3)) - KJK_pol(1d*k,-(i+3),1)

	res = sqrt(-Vinp(hfn,cnr)/[Vinp(hfn,sdx),Vinp(hfn,sdy)])
	if keyword_set(cor) then res = res/sqrt(12*Vinp(hfn,cnr)*product(res))

	return,Cast(res,4,Calctype(k),/fix)
end