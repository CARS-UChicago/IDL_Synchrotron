Function KJK_arg, x, par

;+
; NAME:
;	KJK_ARG
; PURPOSE:
;	Used as kernel in calculation of the KJK function (see there).
; CATEGORY:
;	Mathematical, Integration Kernel.
; CALLING SEQUENCE:
;	Result = KJK_ARG (X, PAR)
; INPUTS:
;    X
;	Numeric, otherwise arbitrary.
;    PAR
;	3 element vector containing [k,/gamma*/theta,/gamma*/psi], or 
;	4 element vector containing the above plus polsub.  Polsub specifies
;	whether pole subtraction should be performed (1) or not (0).  Default
;	setting is 0.
; OPTIONAL INPUT PARAMETERS:
;	None.
; KEYWORD PARAMETERS:
;	None.
; OUTPUTS:
;	Returns the integration kernel for the K.J.Kim function.
; OPTIONAL OUTPUT PARAMETERS:
;	None.
; COMMON BLOCKS:
;	None.
; SIDE EFFECTS:
;	None.
; RESTRICTIONS:
;	None.
; PROCEDURE:
;	Follows the formula in K.J.Kim's paper, NIM A246 (1986) 67.
; MODIFICATION HISTORY:
;	Created 30-MARCH-1994 by Mati Meron.
;-

    tpar = abs([par,0.])
    if tpar(3) and tpar(1) lt tpar(0) then $
    polsub = sqrt(1. - (tpar(1)/tpar(0))^2) else polsub = 0.
    d = (tpar(0)*cos(x) - tpar(1))^2 + (1. + tpar(2)^2)
    sinx = sin(x)

    return, ((d - 2)^2 + 4*tpar(2)^2)*sinx*(sinx - polsub)/d^5
end
